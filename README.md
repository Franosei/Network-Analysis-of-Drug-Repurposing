# A Graph–Bayesian Framework for Drug Repurposing: Integrating Network Centralities, Literature Mining, and Machine Learning

This project identifies drug repurposing opportunities across therapeutic areas by integrating real-world clinical trial data, graph-based centrality metrics, machine learning, and literature-informed Bayesian inference. It uses data from [ClinicalTrials.gov](https://clinicaltrials.gov/) and [PubMed](https://pubmed.ncbi.nlm.nih.gov/) to build and evaluate a graph-based model of drug–disease relationships.

## Objective

To characterise the repurposing potential of drugs by:

- Constructing a bipartite graph of drug–disease relationships from interventional trial data.
- Projecting the graph into drug–drug and disease–disease networks.
- Computing node centrality and structural importance.
- Predicting clinical trial success probabilities using machine learning.
- Integrating semantic literature priors into a Bayesian inference model to rank novel repurposing candidates.

## Project Structure

```
project-root/
├── code/
│   ├── build_training_dataset.py
│   ├── condition_drug_pairs.py
│   ├── data_extraction.py
│   ├── network_builder.py
│   ├── pubmed_utils.py
│   ├── train_phase_success_model.py
│   └── visualize_drug_network.py
├── data/
├── figures/
├── graph/
├── literatures/
├── mesh_data/
│   └── desc2025.xml
├── networks/
├── processed_data/
├── reports/
│   ├── Bayesian_inference_drug_repurposing.pdf
│   └── Network Analysis of Drug Repurposing.pdf
├── .env
├── .gitignore
├── LICENSE
├── note.ipynb
├── README.md
└── requirements.txt
```

## Methodology Overview

1. **Clinical Trial Filtering**  
   - Source: ClinicalTrials.gov v2 API  
   - Filters applied:  
     - Intervention Type: `Drug`  
     - Study Status: `Recruiting`, `Active`, `Completed`  
     - Phase: `Phase 2`, `Phase 3`, `Phase 4`  
   - Output: JSON files per disease area.

2. **Term Normalisation and Mapping**  
   - MeSH 2025 descriptors used for standardisation.  
   - Fuzzy and exact matches resolve inconsistent naming.  
   - Placebos, devices, and unmatched terms excluded.

3. **Network Construction**  
   - **Bipartite Graph**: Nodes = Drugs and Diseases  
   - **Drug–Drug Projection**: Edge weight = shared diseases  
   - Graph metrics computed:  
     - Degree Centrality  
     - Eigenvector Centrality  
     - Betweenness Centrality  
     - Clustering Coefficient  
     - Random Walk (PageRank)

4. **Machine Learning Phase Success Prediction**  
   - Labels based on highest trial phase for each (drug, disease) pair.  
   - Models evaluated:  
     - Logistic Regression  
     - Random Forest (best performance: 77.5%)  
     - XGBoost  
   - Feature importances used as likelihood weights in Bayesian update.

5. **Literature Mining and Semantic Priors**  
   - PubMed queried using drug–disease combinations.  
   - Article conclusions classified as:  
     - Therapeutic  
     - Adverse  
     - Irrelevant  
   - Two priors calculated:  
     - Raw Prior = therapeutic ratio  
     - Penalised Prior = weighted difference of therapeutic and adverse ratios

6. **Bayesian Inference**  
   - Combines penalised prior with network likelihoods.  
   - Outputs posterior distribution over success probability.  
   - Metrics computed:  
     - KL Divergence (information gain)  
     - Posterior Mean Shift (directional change)  
   - Ranks repurposing candidates accordingly.

7. **Visualisation and Reporting**  
   - Network subgraphs of high-centrality drugs  
   - Posterior distribution curves  
   - Centrality heatmaps and ranked candidate tables  
   - All plots saved to `figures/` and `reports/`

## Plots and Outputs Included

- Bipartite Drug–Disease Network GraphML
- Drug–Drug Projected Network with Centrality Scores
- Posterior Probability Distributions
- KL Divergence and Posterior Shift Metrics
- Ranked Table of Candidate Drug–Disease Pairs
- Candidate Literature Abstracts and Evidence Tags
- Final reports in PDF (`reports/`)

## Requirements

- Python 3.8+
- Install dependencies with:

```bash
pip install -r requirements.txt

```

### 2025 MeSH descriptors download:

#### For macOS / Linux / Git Bash:

```bash
mkdir -p mesh_data
wget -P mesh_data/ https://nlmpubs.nlm.nih.gov/projects/mesh/MESH_FILES/xmlmesh/desc2025.xml

```

#### For Windows PowerShell:

```bash
New-Item -ItemType Directory -Force -Path mesh_data
Invoke-WebRequest -Uri https://nlmpubs.nlm.nih.gov/projects/mesh/MESH_FILES/xmlmesh/desc2025.xml -OutFile mesh_data\desc2025.xml

```

- This project is licensed under the MIT License.
© 2025 Francis Osei



