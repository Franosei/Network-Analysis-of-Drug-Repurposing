# _archive/

These scripts are not called by the master runner (`run_full_data_quality_pipeline.py`),
the figure generator (`visualisation/make_all_publication_figures.py`),
the table generator (`reporting/make_all_publication_tables.py`),
or the validation module (`validation/validate_publication_run.py`).

They are preserved here for reference but are not part of the publication pipeline.

## Archived files

| File | Original purpose | Reason archived |
|------|-----------------|-----------------|
| `bayesian_predictor_stability.py` | 10-run stability analysis & plots | Diagnostic tool; outputs go to `supplementary_figures/` |
| `sensitivity_cmax_tau.py` | Grid search for (cmax, tau) prior hyperparameters | Diagnostic tool; not in main pipeline |
| `train_model.py` | Feature weight learning (5 classifiers) | Weights are fixed in DEFAULT_WEIGHTS; not re-run each publication run |
| `visualize_drug_network.py` | Bipartite network visualisation | Not part of evidence-quality manuscript output |
| `weight_validator.py` | Bootstrap/noise/ablation weight validation | Diagnostic tool; supplementary only |
| `section32_candidate_evaluator.py` | Shortlist scoring | Superseded by full ledger and Table4 |

## To restore

Move the relevant file back to `code/` and add a call from the appropriate module.
