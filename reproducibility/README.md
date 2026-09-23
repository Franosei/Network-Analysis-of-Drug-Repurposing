# Reproducibility artefacts

Generated data, registry snapshots, LLM outputs and rendered figures remain
ignored under `outputs/` because they are large and time-dependent. A run is
publication-ready only when its own output directory contains the files listed
in `reference_run_manifest.json`, a dependency snapshot, source acquisition
timestamps and a passing validation report.

The manifest in this directory records the expected artefact contract. It does
not certify an old run: the existing 23 September 2026 output predates the
current safety-status and positive-unlabelled validation controls and must be
regenerated before publication.

For an archival release, copy the complete validated run directory to Zenodo,
a GitHub Release or an equivalent repository, then record its checksums in the
manifest. Do not treat a current ClinicalTrials.gov record filtered by date as
a historical snapshot.
