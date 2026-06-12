"""Metric definitions used by the publication-quality audit pipeline."""

from __future__ import annotations


LEDGER_COLUMNS = [
    # Identity / provenance
    "drug_original",
    "disease_original",
    "drug_cleaned",
    "disease_cleaned",
    "drug_mapped",
    "disease_mapped",
    "drug_mesh_id",
    "disease_mesh_id",
    # Retained short-form aliases used in downstream code
    "drug",
    "disease",
    # Terminology mapping audit
    "drug_mapping_method",
    "disease_mapping_method",
    "drug_mapping_score",
    "disease_mapping_score",
    "drug_token_jaccard_score",
    "disease_token_jaccard_score",
    "drug_mapping_status",
    "disease_mapping_status",
    # Clinical trial evidence
    "trial_count",
    "phase_distribution",
    "status_distribution",
    # Literature evidence
    "articles_retrieved",
    "usable_articles",
    "therapeutic_count",
    "adverse_count",
    "irrelevant_count",
    "therapeutic_ratio",
    "adverse_burden",
    "irrelevant_noise_rate",
    "literature_completeness_score",
    # Safety overlap
    "safety_overlap_gamma",
    "safety_penalty",
    # Graph / structural evidence
    "graph_distance",
    "random_walk_score",
    "katz_similarity",
    "preferential_attachment",
    "structural_likelihood",
    "structural_consistency_score",
    # Bayesian inference
    "prior_mean",
    "posterior_mean",
    "credible_interval_lower",
    "credible_interval_upper",
    "credible_interval_width",
    "kl_divergence",
    "mean_shift",
    # Composite scores and flags
    "evidence_readiness_score",
    "uncertainty_level",
    "quality_flag",
    "coverage_tier",
    "final_interpretation",
]

# Shorter alias list used for backward-compatible validation — subset that must
# always be present in a ledger saved to disk.
REQUIRED_LEDGER_COLUMNS = {
    "drug",
    "disease",
    "drug_mapping_method",
    "disease_mapping_method",
    "drug_mapping_score",
    "disease_mapping_score",
    "drug_mapping_status",
    "disease_mapping_status",
    "trial_count",
    "articles_retrieved",
    "therapeutic_count",
    "adverse_count",
    "irrelevant_count",
    "safety_overlap_gamma",
    "posterior_mean",
    "credible_interval_width",
    "kl_divergence",
    "evidence_readiness_score",
    "quality_flag",
    "coverage_tier",
}


METRIC_DEFINITIONS = {
    "evidence_readiness_score": (
        "Composite 0-100 score for evidence quality and audit readiness, not clinical efficacy. "
        "Weights: entity_mapping=0.15, literature_completeness=0.15, semantic_relevance=0.15, "
        "adverse_cleanliness=0.10, safety_alignment=0.15, structural_consistency=0.15, "
        "posterior_certainty=0.15."
    ),
    "coverage_tier": (
        "Completeness of local evidence coverage for a drug-disease pair. "
        "Values: full_bayesian_audit, bayesian_without_safety, literature_and_graph, "
        "graph_only, literature_only, matched_pairs_only."
    ),
    "therapeutic_ratio": "Therapeutic-classified articles divided by all usable classified articles.",
    "adverse_burden": "Adverse-classified articles divided by all usable classified articles.",
    "irrelevant_noise_rate": "Irrelevant classified articles divided by all classified articles.",
    "literature_completeness_score": (
        "0-1 score reflecting article retrieval completeness: penalises zero-article pairs, "
        "rewards pairs where usable articles are a high fraction of retrieved articles."
    ),
    "safety_overlap_gamma": (
        "Semantic overlap score between adverse-event terms and the target disease "
        "symptoms/phenotype. Range [0, 1]; higher values indicate stronger adverse overlap."
    ),
    "safety_penalty": (
        "Multiplicative penalty applied to the prior mean due to adverse-event/disease overlap. "
        "p_penalised = p_raw × (1 − penalty_scale × gamma)."
    ),
    "structural_consistency_score": (
        "Normalised composite of graph-feature signals (graph_distance, random_walk_score, "
        "katz_similarity, preferential_attachment, structural_likelihood). Range [0, 1]."
    ),
    "credible_interval_width": (
        "Width of the posterior 95% credible interval (Beta distribution). "
        "Narrow intervals indicate well-constrained posteriors; wide intervals indicate uncertainty."
    ),
    "kl_divergence": (
        "Kullback-Leibler divergence from prior to posterior Beta distributions. "
        "Measures information gain; low KL means the data barely updated prior beliefs."
    ),
    "mean_shift": (
        "Posterior mean minus prior mean. Positive values indicate evidence pushed belief "
        "toward repurposing; negative values indicate evidence against."
    ),
    "uncertainty_level": (
        "Categorical label derived from credible_interval_width: "
        "low (<0.15), moderate (0.15-0.35), high (>0.35)."
    ),
    "drug_mapping_method": (
        "How the drug term was normalised to MeSH. "
        "Values: exact_match, fuzzy_high_confidence (score≥0.80), "
        "fuzzy_low_confidence (score<0.80), unmapped, mapped_legacy_no_provenance."
    ),
    "disease_mapping_method": (
        "How the disease term was normalised to MeSH. "
        "Values: exact_match, fuzzy_high_confidence (score≥0.80), "
        "fuzzy_low_confidence (score<0.80), unmapped, mapped_legacy_no_provenance."
    ),
    "final_interpretation": (
        "Human-readable audit summary: qualitative assessment of the drug-disease pair's "
        "evidence-readiness status, flagging any unresolved concerns."
    ),
}
