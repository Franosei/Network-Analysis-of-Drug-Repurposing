# Evidence-Readiness Network Analysis for Drug Repurposing

This repository implements an auditable drug-repurposing pipeline that combines ClinicalTrials.gov trial structure, MeSH terminology standardisation, interpretable graph features, PubMed/PMC semantic evidence, openFDA safety overlap, and Bayesian uncertainty quantification.

The output is not a clinical efficacy claim. It is an evidence-readiness ledger: a technical audit of whether a drug-disease pair has enough clean, concordant, structurally plausible, and uncertainty-bounded evidence to justify expert review.

## Pipeline Diagram

![Evidence-readiness assessment pipeline for drug repurposing](pipeline/pipeline_drug_repurposing.png)

The diagram summarises the end-to-end publication workflow:

1. ClinicalTrials.gov drug-intervention studies are filtered, deduplicated, and expanded into drug-condition rows.
2. Drug and disease terms are normalised against MeSH descriptors and entry terms.
3. Pairs are assigned an evidence coverage tier: matched only, graph only, literature only, literature plus graph, or full Bayesian audit.
4. Literature and safety evidence are retrieved and semantically classified.
5. A bipartite drug-disease graph produces structural likelihood features.
6. Semantic priors and graph likelihoods are fused into posterior Beta distributions.
7. Composite evidence-readiness scores and rule-based quality categories are written to a per-pair ledger.

## Core Entry Points

| File | Purpose |
| --- | --- |
| `run_full_data_quality_pipeline.py` | Main orchestration script. Creates a dated publication run folder with numbered audit files, ledgers, tables, figures, logs, and validation output. |
| `code/data_extraction.py` | ClinicalTrials.gov v2 API fetcher with paging, retry/backoff, status/phase/intervention filters, cutoff-year filtering, and global NCT de-duplication. |
| `code/condition_drug_pairs.py` | Extracts condition-intervention pairs and maps terms to MeSH using exact, fuzzy, and token-guided matching. |
| `code/scalable_graph_builder.py` | Canonical leakage-controlled graph builder. Target edges are withheld before feature calculation and sampled non-edges are labelled as unlabelled, not factual negatives. |
| `code/network_builder.py` | Legacy exploratory graph builder retained for historical comparison only. It is not used by the publication pipeline because it can leak observed edges into predictors. |
| `code/train_model.py` | Legacy graph-weight experiments retained for historical comparison only. Its row-level known/unknown split is not a publication metric. |
| `code/pubmed_utils.py` | Retrieves PubMed records, optionally augments with PMC text, classifies abstracts with an LLM, and constructs literature priors. |
| `code/side_effect_updater.py` | Retrieves openFDA adverse-event terms and applies an LLM-scored safety-overlap penalty to the prior. |
| `code/bayesian_predictor.py` | Fuses literature/safety priors with graph likelihoods using Beta pseudo-count updates and writes run logs/plots. |
| `code/evidence_quality_report.py` | Offline report generator that builds pair-level, source-level, and summary quality tables from existing artifacts. |
| `data_quality/*.py` | Shared metric definitions, composite scoring weights, coverage-tier rules, and quality-flag rules. |
| `reporting/make_all_publication_tables.py` | Generates manuscript and supplementary CSV tables. |
| `visualisation/make_all_publication_figures.py` | Generates publication and supplementary figures. |
| `validation/validate_publication_run.py` | Validates a publication run folder for required structure, ledger columns, panel coverage, tables, figures, and reproducibility snapshots. |

## Repository Layout

```text
.
|-- code/                         # Core extraction, mapping, graph, literature, Bayesian, and utility modules
|-- config/                       # Case-study panels used by publication runs
|-- data_quality/                 # Ledger schema, scoring weights, and rule-based flags
|-- outputs/                      # Ignored local publication/data-quality rerun folders
|   `-- <run>/                     # Generated locally or restored from an archive
|-- pipeline/                     # Pipeline documentation figure
|-- reporting/                    # Publication table generation
|-- validation/                   # Publication run validation
|-- visualisation/                # Publication figure generation
|-- run_full_data_quality_pipeline.py
|-- requirements.txt
`-- README.md
```

## Methodology

### 1. Clinical Trial Extraction

`ClinicalTrialFetcher` queries the ClinicalTrials.gov v2 API by therapeutic area and filters studies locally:

- `overallStatus` in `RECRUITING`, `ACTIVE_NOT_RECRUITING`, `ENROLLING_BY_INVITATION`, `COMPLETED`, `WITHDRAWN`, or `TERMINATED`
- phase intersects `PHASE2`, `PHASE3`, or `PHASE4`
- at least one intervention has type `DRUG`
- reference year is at or before `--cutoff_year` (default `2020`)
- duplicate `NCTId` values are removed globally across therapeutic areas

The publication pipeline writes `01_clinical_trial_extraction_audit.csv` with source-level counts and filtering diagnostics.

### 2. Terminology Standardisation

`ConditionDrugPairBuilder` loads `mesh_data/desc2026.xml` and builds normalised lookup maps from descriptor names and entry terms:

- disease terms are identified from MeSH tree numbers beginning with `C`
- drug terms are identified from MeSH tree numbers beginning with `D`
- incoming terms retain clinically meaningful modifiers such as stage, severity, phenotype, response state, and positive/negative status; raw terms, normalised terms, match method, confidence, MeSH ID, and ambiguity metadata are retained
- matching proceeds through exact match, high-confidence fuzzy match, low-confidence fuzzy match, and token-guided Jaccard scoring

In a publication run, mapped pairs are written under `outputs/<run>/processed_data/`; unresolved terms and failure reasons are written to `outputs/<run>/processed_data/unmatched_pairs.json`. The mapping audit is mirrored into `02_terminology_mapping_audit.csv`.

### 3. Graph Construction and Structural Evidence

`scalable_graph_builder.py` is the canonical leakage-controlled implementation. The historical `InterpretableGraphFeatureBuilder` is retained only for exploratory comparison and is not part of the publication pathway.

- drug nodes connect to disease nodes through observed trial pairs
- held-out observed pairs receive `Label = 1` while their target edge is removed before feature calculation
- sampled non-edges are marked `LabelSemantics = unlabelled_non_edge`; absence from the registry is not treated as a factual negative
- the graph is written as GraphML for downstream inspection

The reported ROC-AUC, PR-AUC, precision-at-k and recall-at-k are for a
transductive edge-reconstruction sensitivity analysis. Drugs and diseases may
occur on both sides of that split. They must not be described as prospective
clinical prediction, drug-discovery performance, temporal generalisation or
efficacy validation. Drug-disjoint, disease-disjoint and temporal evaluations
remain separate future validation designs rather than hidden assumptions.

For each pair, the graph stage computes five interpretable features:

| Feature | Interpretation |
| --- | --- |
| `GraphDistanceToIndication` | `1/3` when at least one alternative length-three bipartite path exists after withholding the target edge, otherwise `0`. |
| `RandomWalkScore` | Degree-normalised contribution of the observed length-three alternative paths. |
| `StructuralLikelihood` | `(1 + drug centrality) * (1 + disease centrality)` using degree centrality on the withheld graph. |
| `PreferentialAttachment` | Product of the two degree-centrality values. |
| `KatzSimilarity` | `0.005^3 *` the number of alternative length-three paths. |

Outputs inside a publication run:

- `outputs/<run>/graph/bipartite_target_edge_withheld.graphml`
- `outputs/<run>/graph/graph_features_known.csv`
- `outputs/<run>/graph/graph_features_unknown.csv` (requested target pairs that were not observed)
- `03_graph_construction_audit.csv` inside publication runs

### 4. Literature Prior and Safety Adjustment

`LLMClassifier.build_semantic_prior()` retrieves PubMed records for each drug-disease pair, deduplicates PMIDs, classifies usable records, and saves classified evidence JSON.

Each record is labelled as one mutually exclusive evidence direction:

- `benefit`
- `null`
- `harm`
- `conflicting`
- `irrelevant`

The legacy `therapeutic` and `adverse` totals remain as derived compatibility fields. Study design, DOI, journal, PMID and publication-type metadata are retained, and repeated evidence units are deduplicated before concentration is assigned.

Given `T` benefit records, `A` records in the null, harm or conflicting classes, and `K` classified records:

```text
p_raw        = T / K
p_penalised  = max((T - 2A) / K, 0)
```

The separate concentration input `M` is the number of deduplicated evidence units and is used only for the pseudo-count strength `c(M)`. Thus repeated records do not narrow the interval merely by increasing the number of retrieved rows.

The adverse-evidence coefficient `2` is a pre-specified conservative heuristic: one adverse-classified article is allowed to offset two therapeutic-classified articles before clamping at zero. This intentionally penalises safety-conflicted literature more strongly than mixed neutral evidence.

`SideEffectUpdater` then retrieves openFDA adverse-event terms and asks the LLM whether those terms semantically overlap with the target disease. If a relation exists:

```text
p_final = p_penalised * (1 - penalty_scale * gamma)
```

where `gamma` is the safety-overlap confidence in `[0, 1]` and `penalty_scale` defaults to `0.5`. Missing safety data are represented as `gamma = null` and `safety_data_status = API_ERROR`, `NO_RECORDS`, `LLM_ERROR`, or `API_KEY_MISSING`; they are never interpreted as evidence of no safety concern. When the openFDA response reaches the configured report limit, the status is `COMPLETE_LIMITED` to make the sampling boundary explicit.

The safety penalty scale `0.5` is also a pre-specified conservative heuristic. It caps the safety-overlap reduction at 50% when `gamma = 1`, avoiding complete elimination of a signal while still down-weighting candidates whose adverse-event profile overlaps with the target phenotype.

Publication runs write:

- `04_literature_retrieval_audit.csv`
- `05_semantic_classification_audit.csv`
- `06_safety_overlap_audit.csv`

### 5. Bayesian Fusion

The Bayesian stage converts the safety-adjusted literature prior and graph likelihood into an Evidence Posterior Index using Beta pseudo-counts. The interval is a pseudo-Beta uncertainty interval. It is not a calibrated 95% clinical efficacy interval.

Evidence-scaled prior concentration:

```text
c(M) = cmax * (1 - exp(-M / tau))
```

Default values:

```text
cmax               = 200
tau                = 25
likelihood_strength = 50
likelihood_intercept = 0
```

Prior:

```text
alpha_prior = 1 + c(M) * p_final
beta_prior  = 1 + c(M) * (1 - p_final)
```

Graph likelihood:

```text
score  = intercept + sum(weight_i * feature_i)
p_like = sigmoid(score)
```

Default graph-feature weights in `run_full_data_quality_pipeline.py`:

```text
GraphDistanceToIndication =  0.1148
RandomWalkScore           =  0.2470
StructuralLikelihood      = -0.1154
PreferentialAttachment    = -0.0410
KatzSimilarity            =  1.6515
```

Likelihood pseudo-counts:

```text
alpha_like = 1 + likelihood_strength * p_like
beta_like  = 1 + likelihood_strength * (1 - p_like)
```

Posterior fusion subtracts the unit baselines and adds only evidence mass:

```text
alpha_post = alpha_prior + (alpha_like - 1)
beta_post  = beta_prior  + (beta_like  - 1)
posterior_mean = alpha_post / (alpha_post + beta_post)
```

Diagnostics include posterior variance, 95% credible interval, credible-interval width, KL divergence from prior to posterior, and posterior mean shift.

The Bayesian constants above are fixed defaults used by the pipeline and are written to `run_config.json` for each publication run.

### 6. Composite Evidence-Readiness Score

The final 0-100 score is a weighted audit-readiness index, not a treatment recommendation:

```text
score = 100 * (
    0.15 * entity_mapping
  + 0.15 * literature_completeness
  + 0.15 * semantic_relevance
  + 0.10 * adverse_cleanliness
  + 0.15 * safety_alignment
  + 0.15 * structural_consistency
  + 0.15 * posterior_certainty
)
```

Rule-based quality categories include:

- `High evidence quality`
- `Moderate evidence quality`
- `Low evidence quality`
- `Terminology uncertainty`
- `Insufficient evidence`
- `Sparse literature; structurally plausible`
- `Literature noise dominated`
- `Literature-conflicted`
- `Safety-concerning`
- `Safety-conflicted evidence`

Coverage tiers are standardised as:

- `full_bayesian_audit`
- `bayesian_without_safety`
- `literature_and_graph`
- `graph_only`
- `literature_only`
- `matched_pairs_only`
- `insufficient_coverage`

## Main Outputs

The retained full run is `outputs/20260610_bayesian/`. A publication run folder has this structure:

```text
outputs/<run>/
|-- audit_files/
|   |-- 01_clinical_trial_extraction_audit.csv
|   |-- 02_terminology_mapping_audit.csv
|   |-- 03_graph_construction_audit.csv
|   |-- 04_literature_retrieval_audit.csv
|   |-- 05_semantic_classification_audit.csv
|   |-- 06_safety_overlap_audit.csv
|   |-- 07_bayesian_uncertainty_audit.csv
|   |-- 08_full_evidence_quality_ledger.csv
|   |-- 09_quality_category_counts.csv
|   |-- 10_summary_dashboard.csv
|   |-- 11_disease_drug_pair_validation.csv
|   `-- 12_case_study_validation.csv
|-- graph/
|-- ledgers/
|   `-- full_evidence_quality_ledger.csv
|-- logs/
|   |-- run_log.txt
|   `-- requirements_snapshot.txt
|-- manuscript_figures/
|-- manuscript_tables/
|-- supplementary_figures/
|-- supplementary_tables/
|-- run_config.json
|-- validation_checks.csv
`-- validation_report.json
```

The central artifact is:

```text
outputs/<run>/ledgers/full_evidence_quality_ledger.csv
```

Key ledger columns include entity provenance, MeSH IDs, mapping scores, trial counts, literature counts, safety overlap, graph features, prior/posterior summaries, credible intervals, KL divergence, composite readiness score, coverage tier, and final interpretation.

## Installation

Python 3.9 to 3.12 is the supported range for the pinned dependency set in
`requirements.txt`. Record the exact interpreter and installed versions in each
publication run's `logs/requirements_snapshot.txt`.

```powershell
py -m venv .venv
.\.venv\Scripts\Activate.ps1
py -m pip install --upgrade pip
py -m pip install -r requirements.txt
```

The Bayesian literature and safety stages require an OpenAI API key. Put it in `.env`:

```text
OPENAI_API_KEY=your_key_here
```

ClinicalTrials.gov, PubMed, PMC, and openFDA are queried over public HTTP APIs. Full reruns therefore require network access.

`code/train_model.py` imports `xgboost`; install it separately if you want to rerun the optional model-comparison training stage.

## MeSH Descriptor File

The mapper expects:

```text
mesh_data/desc2026.xml
```

Download it with PowerShell:

```powershell
New-Item -ItemType Directory -Force -Path mesh_data
Invoke-WebRequest -Uri https://nlmpubs.nlm.nih.gov/projects/mesh/MESH_FILES/xmlmesh/desc2026.xml -OutFile mesh_data\desc2026.xml
```

Or let the full pipeline download it:

```powershell
py run_full_data_quality_pipeline.py --download_mesh true
```

## Running the Pipeline

### Full publication rerun

This refreshes ClinicalTrials.gov extraction, MeSH mapping, graph construction, PubMed/PMC evidence, safety overlap, Bayesian scoring, tables, figures, and validation.

```powershell
py run_full_data_quality_pipeline.py `
  --full_registry true `
  --download_mesh true `
  --mesh_exact_only true `
  --pubmed_max_articles 0 `
  --pubmed_filter_level exact `
  --pubmed_years_back 0 `
  --run_bayesian true `
  --refresh_literature true `
  --refresh_safety true `
  --output_dir outputs/publication_run_YYYYMMDD
```

### Faster structural/data-quality rerun

This refreshes clinical extraction, mapping, and graph artifacts but reuses existing literature/safety/Bayesian artifacts.

```powershell
py run_full_data_quality_pipeline.py `
  --run_bayesian false `
  --refresh_literature false `
  --refresh_safety false `
  --output_dir outputs/publication_structural_YYYYMMDD
```

### Offline report from existing artifacts

```powershell
py code\evidence_quality_report.py `
  --matched outputs\20260610_bayesian\processed_data\condition_drug_pairs.json `
  --unmatched outputs\20260610_bayesian\processed_data\unmatched_pairs.json `
  --known-graph outputs\20260610_bayesian\graph\graph_features_known.csv `
  --unknown-graph outputs\20260610_bayesian\graph\graph_features_unknown.csv `
  --runs-dir outputs\20260610_bayesian\runs `
  --literature-dir outputs\20260610_bayesian\literatures `
  --output-dir outputs\20260610_bayesian\manuscript_tables
```

### Generate tables or figures manually

```powershell
py reporting\make_all_publication_tables.py `
  --ledger_path outputs\publication_run_YYYYMMDD\ledgers\full_evidence_quality_ledger.csv `
  --audit_dir outputs\publication_run_YYYYMMDD\audit_files `
  --output_dir outputs\publication_run_YYYYMMDD\manuscript_tables `
  --supp_dir outputs\publication_run_YYYYMMDD\supplementary_tables `
  --panel_csv config\case_study_panel_publication.csv
```

```powershell
py visualisation\make_all_publication_figures.py `
  --ledger_path outputs\publication_run_YYYYMMDD\ledgers\full_evidence_quality_ledger.csv `
  --audit_dir outputs\publication_run_YYYYMMDD\audit_files `
  --runs_dir outputs\publication_run_YYYYMMDD\runs `
  --output_dir outputs\publication_run_YYYYMMDD\manuscript_figures `
  --supp_dir outputs\publication_run_YYYYMMDD\supplementary_figures `
  --panel_csv config\case_study_panel_publication.csv
```

### Generate supplementary sensitivity tables

This offline analysis varies `cmax`, `tau`, likelihood strength `lambda`, the adverse-evidence weight, and the safety coefficient by 50% below and above their baseline values. It holds ClinicalTrials.gov, MeSH, PubMed/PMC, safety-gamma, and graph artifacts fixed.

```powershell
py reporting\make_sensitivity_supplement.py `
  --run_dir outputs\20260610_bayesian
```

Outputs:

- `outputs/20260610_bayesian/supplementary_tables/SuppTable_sensitivity_parameter_summary.csv`
- `outputs/20260610_bayesian/supplementary_tables/SuppTable_sensitivity_pair_level.csv`
- `outputs/20260610_bayesian/supplementary_tables/SuppText_sensitivity_analysis.md`

### Validate a publication run

```powershell
py validation\validate_publication_run.py `
  --output_dir outputs\publication_run_YYYYMMDD `
  --panel_csv config\case_study_panel_publication.csv
```

Use `--strict` to return a non-zero exit code when validation errors are present.

## Full HCQ to COVID-19 rerun: methodology and results

The complete hydroxychloroquine and COVID-19 case study was rerun on 23 September 2026. This rerun replaced the former limited literature retrieval with an open retrieval setting. In the current pipeline, `--pubmed_max_articles 0` means that every result returned by the declared PubMed query is retrieved. There is no fixed article-count target.

### Data cut and acquisition (reference run only)

The values below describe the archived reference run, not a promise that an
online source will return the same records later. Every new run must record its
own UTC acquisition time, query, filters, source status and artifact paths in
`run_config.json` and the stage audit files.

| Source | Acquisition time | Data used |
|---|---|---|
| ClinicalTrials.gov v2 | 23 September 2026, 13:27 to 13:37 UTC | Current registry snapshot, all studies |
| MeSH | `desc2026.xml` | Descriptor names, entry terms, and MeSH identifiers |
| PubMed | 23 September 2026, 13:49 UTC | Exact title and abstract query, with no date restriction or article cap |
| openFDA FAERS | During the Bayesian run, approximately 15:50 to 16:02 UTC | Hydroxychloroquine adverse-event terms |

The registry and PubMed results are time-dependent. The reported values describe these snapshots and should not be presented as permanent database totals.

### 1. ClinicalTrials.gov extraction

The ClinicalTrials.gov v2 API was queried without therapeutic-area sampling or a year cut-off. The selected fields were NCT identifier, title, status, start and posting dates, conditions, interventions, intervention type, phase, and study type.

Studies were retained for graph construction when they were interventional, contained at least one intervention explicitly typed as `DRUG`, and contained at least one condition. Placebo-only interventions were excluded. Drug and condition fields were crossed within each study, retaining the NCT identifier and provenance.

The full extraction contained:

- 604,146 registry studies;
- 201,229 interventional studies with drug interventions;
- 680,694 raw drug-condition study rows;
- 61 strict source-term hydroxychloroquine and COVID-19 trials.

The strict target trials were first posted from 7 February 2020 onwards. This is distinct from the MeSH-canonical count below, because MeSH combines equivalent source terms.

### 2. MeSH harmonisation

Terms were mapped to MeSH descriptors and entry terms using exact normalised matching. The normalisation is case-insensitive and removes punctuation and spacing variation, but it does not apply fuzzy matching or unsupported synonyms in the exact full-registry run.

The target mappings were:

```text
Hydroxychloroquine  -> D006886
COVID-19             -> D000086382
```

The mapped target pair contained 94 distinct NCT trials after equivalent source terms were harmonised. The mapped registry contained 139,837 pair rows, 36,598 unique canonical drug-condition pairs, and 167,604 unresolved entity records retained for audit.

The distinction is therefore:

```text
61  strict source-term HCQ-COVID trials
94  MeSH-canonical HCQ-COVID trials
```

### 3. Full graph and weight update

A bipartite graph was built from all observed MeSH-canonical drug-condition pairs. The direct hydroxychloroquine and COVID-19 edge was withheld before target scoring so that the target trial records could not directly inflate its graph evidence.

Graph weights were refitted using balanced logistic regression. Twenty per cent of registered edges were held out for training examples, an equal number of non-edges were sampled, and the target pair was excluded from weight fitting. The held-out test ROC-AUC was 0.9339.

The updated raw feature weights were:

```text
GraphDistanceToIndication    0.8364
RandomWalkScore            706.8406
StructuralLikelihood         79.8871
PreferentialAttachment   -10507.9407
KatzSimilarity            91836.2090
Intercept                   -82.3092
```

For the withheld target pair, the graph features were:

```text
Graph distance                 0.3333
Random-walk score              0.004953
Structural likelihood          1.0663
Preferential attachment        0.000997
Katz similarity                0.000283
Graph likelihood probability    0.999999
```

The graph likelihood represents registry co-registration and network topology. It is not a clinical efficacy estimate.

### 4. Open PubMed retrieval and classification

The declared query was:

```text
"hydroxychloroquine"[Title/Abstract] AND "COVID-19"[Title/Abstract]
```

The query had no date restriction and no article cap. PubMed returned 3,391 identifiers, and all 3,391 records were parsed. Exact phrase verification retained 3,108 records containing both terms in the title or abstract. The remaining 283 records were retained in the retrieval audit but excluded from the exact-pair evidence set.

The 3,108 exact-pair records were classified with the configured `gpt-4o-mini` semantic classifier using five mutually exclusive classes: benefit, null, harm, conflicting and irrelevant. For the compact reference summary, null, harm and conflicting are combined into the derived adverse/conflicting total:

```text
Benefit                  893
Null, harm or conflicting 1271
Irrelevant               944
Total                   3108
```

### 5. Safety adjustment

Current FAERS terms for hydroxychloroquine were compared with COVID-19 clinical features using the configured semantic safety check. The run returned eight matching terms and a safety-overlap value of:

```text
gamma = 0.75
```

The pre-specified safety penalty scale was 0.5. This safety result is a reported-signal overlap and does not establish causality.

### 6. Bayesian fusion

The literature prior was calculated as:

```text
p_raw = T / M
      = 893 / 3108
      = 0.2873
```

The adverse-evidence penalty was:

```text
p_penalised = max((T - 2A) / M, 0)
             = 0.000001
```

After the safety adjustment:

```text
p_final = 0.000001
```

The Bayesian constants were `cmax = 200`, `tau = 25`, likelihood strength `50`, and likelihood intercept `0`. With 3,108 classified records, the evidence concentration reached `c(M) = 200`.

The resulting distributions were:

```text
Prior:     Beta(1.0002, 200.9998)
Likelihood Beta(50.99995, 1.00005)
Posterior: Beta(51.00015, 200.99985)
```

### Final result

| Measure | Result |
|---|---:|
| Posterior mean | 0.2024 |
| 95% credible interval | 0.1552 to 0.2540 |
| Posterior variance | 0.000638 |
| Credible-interval width | 0.0989 |
| KL divergence | 29.6651 |
| Evidence-readiness score | 74.872 |
| Quality classification | Safety-conflicted evidence |
| Coverage tier | Full Bayesian audit |

The Evidence Posterior Index is a model-based audit signal under the stated pseudo-count construction. Its interval is not a calibrated probability interval for clinical efficacy. The high graph likelihood is outweighed by the adverse or conflicting literature count and the safety-overlap penalty.

### Reproducibility artefacts

The tracked artefact contract is defined in
`reproducibility/reference_run_manifest.json`. Generated outputs are ignored
because registry, PubMed, FAERS and LLM results are time-dependent. A clean
rerun must pass the current validator before its checksums are added to that
manifest and the run is archived.

- `outputs/20260923_exact_pair_hcq_covid_full_registry/run_manifest.json`: registry and PubMed acquisition manifest.
- `outputs/20260923_hcq_covid_full_pipeline/processed_data/condition_drug_pairs.json`: MeSH-mapped trial pairs.
- `outputs/20260923_hcq_covid_full_pipeline/graph/updated_graph_weights.json`: refitted graph weights.
- `outputs/20260923_hcq_covid_full_pipeline/graph/graph_features_known.csv`: withheld-edge target graph features.
- `outputs/20260923_hcq_covid_full_pipeline/runs/run_hydroxychloroquine_covid-19_20260923_155027.json`: Bayesian components and posterior parameters.
- `outputs/20260923_hcq_covid_full_pipeline/ledgers/full_evidence_quality_ledger.csv`: final evidence ledger.

To reproduce the open literature setting in a new full-registry run, use:

```powershell
py run_full_data_quality_pipeline.py `
  --full_registry true `
  --full_registry_db outputs\20260923_exact_pair_hcq_covid_full_registry\clinicaltrials_exact_graph.sqlite `
  --resume_full_registry true `
  --mesh_exact_only true `
  --pubmed_max_articles 0 `
  --pubmed_filter_level exact `
  --pubmed_years_back 0 `
  --refresh_literature true `
  --refresh_safety true `
  --run_bayesian true `
  --output_dir outputs\full_hcq_covid_rerun
```

## Case-Study Panel

The publication panel is defined in:

```text
config/case_study_panel_publication.csv
```

It includes successful, failed, conflicted, emerging, noisy, and sparse-evidence examples. The panel is used to test whether the evidence-readiness framework distinguishes high-quality repurposing evidence from literature volume, structural plausibility, terminology uncertainty, and safety conflict.

## Validation Philosophy

The validation layer checks reproducibility artifacts and audit completeness rather than model accuracy alone. It verifies:

- required output directories
- `run_config.json`
- run logs and dependency snapshots
- required ledger columns
- duplicate drug-disease rows
- case-study panel coverage
- manuscript tables and figures
- row-count consistency

Warnings are expected when a run intentionally reuses legacy artifacts, lacks fresh literature, lacks safety overlap, or does not include every case-study pair in the active ledger. Errors should be resolved before treating a run as publication-ready.

## Modelling Limitations

The adverse-evidence coefficient `2`, safety penalty scale `0.5`, prior concentration parameters `cmax=200` and `tau=25`, and graph likelihood strength `50` are transparent, pre-specified modelling constants rather than fitted causal parameters. They make the current framework conservative for noisy or safety-conflicted evidence, but future work should quantify sensitivity to alternative adverse-penalty, safety-penalty, prior-concentration, and likelihood-strength settings.

The semantic labels are generated by an LLM under the recorded model name and
`five-way-evidence-v1` schema. The repository does not claim human-labelled
classification accuracy, sensitivity, specificity or posterior calibration;
those require an independent, frozen validation set. Registry, literature and
network sources can also describe the same underlying studies, so the fused
pseudo-count construction should be read as an auditable aggregation, not as
independent likelihood evidence.

## License

MIT License. See `LICENSE`.

Copyright 2025 Francis Osei.
