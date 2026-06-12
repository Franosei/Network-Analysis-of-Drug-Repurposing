# Data Quality and Evidence Quality Report

This report was generated offline from existing pipeline artifacts. It measures evidence readiness and data quality, not clinical efficacy.

## Outputs

- `pair_level_evidence_quality.csv`: one row per audited drug-disease run log.
- `source_level_data_quality_audit.csv`: ClinicalTrials.gov, MeSH, PubMed/PMC, safety, graph, and posterior quality metrics.
- `summary_dashboard.csv`: high-level counts, rates, and quality category totals.
- `clinical_trial_data_readiness.csv`: extraction/normalisation readiness metrics.
- `terminology_standardisation_quality.csv`: mapping-method proportions and nominal confidence.

## Clinical-Trial Readiness

- Drug-condition rows audited: 31235
- Matched rows: 10178
- Unmatched condition rate: 0.203842
- Unresolved drug rate: 0.590719
- Final graph-ready unique non-placebo pairs: 5026

## Pair-Level Evidence

- Audited pair count: 56
- Average evidence readiness score: 62.209
- Average posterior uncertainty width: 0.199359

Top evidence-readiness pairs:
- Hydroxychloroquine -> Melanoma: 81.447 (High evidence quality)
- Levodopa -> cardiovascular diseases: 77.008 (Literature-conflicted)
- metformin -> hypersensitivity: 76.750 (High evidence quality)
- metformin -> pulmonary embolism: 76.185 (High evidence quality)
- dexamethasone -> COVID-19: 75.964 (Safety-concerning)
- metformin -> mitochondrial diseases: 75.556 (Literature-conflicted)
- metformin -> mania: 75.291 (High evidence quality)
- metformin -> foot ulcer: 74.392 (Moderate evidence quality)
- metformin -> ascites: 72.757 (Moderate evidence quality)
- metformin -> maple syrup urine disease: 71.915 (Insufficient evidence)

## Interpretation Notes

- Evidence readiness scores combine mapping reliability, literature completeness, semantic relevance, adverse burden, safety overlap, structural consistency, and posterior uncertainty.
- Successful MeSH mappings in the current legacy processed pair file do not include exact/fuzzy provenance, so they are labelled `mapped_legacy_no_provenance`.
- Safety overlap terms are included when present in run logs; older logs may only contain gamma.

## Dashboard Snapshot

- report_generated_at: 2026-06-10T11:40:58
- pair_count: 56
- clinical_drug_condition_rows: 31235
- final_graph_ready_pair_count: 5026
- mapping_success_rate: 0.561005
- unresolved_entity_rate: 0.438995
- usable_literature_rate: 1.0
- average_irrelevant_retrieval_rate: 0.34069
- average_safety_overlap_score: 0.35625
- average_posterior_uncertainty: 0.199359
- average_evidence_readiness_score: 62.209
- median_evidence_readiness_score: 61.232
- quality_category_count:High evidence quality: 4
- quality_category_count:Insufficient evidence: 16
- quality_category_count:Literature noise dominated: 1
- quality_category_count:Literature-conflicted: 4
- quality_category_count:Moderate evidence quality: 6
- quality_category_count:Safety-concerning: 4
- quality_category_count:Safety-conflicted evidence: 1
- quality_category_count:Sparse literature; structurally plausible: 3
- quality_category_count:Terminology uncertainty: 17
- evidence_domain_support_count:literature_and_network: 3
- evidence_domain_support_count:literature_only: 21
- evidence_domain_support_count:neither: 29
- evidence_domain_support_count:network_only: 3
- legacy_mapping_provenance_note: Current matched pair artifacts do not store exact vs fuzzy match provenance for successful mappings.
