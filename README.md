# Network Analysis of Drug Repurposing Across Therapeutic Areas

This project leverages clinical trial data from [ClinicalTrials.gov](https://clinicaltrials.gov/) to identify potential drug repurposing opportunities across therapeutic areas using **network science techniques**. The analysis focuses exclusively on active pharmaceutical interventions (excluding devices and procedures), and considers trials in Phase I–IV with valid recruitment statuses.

## Objective

To identify and characterise repurposing opportunities by analysing the co-occurrence of drug interventions across multiple disease categories, using centrality and community detection in a drug–disease network constructed from real-world clinical trials.

## Project Structure

<pre> <code> ## Project Structure ``` ├── code/ │ ├── data_extraction.py # Downloads and filters trials from ClinicalTrials.gov API │ ├── build_graphs.py # Builds bipartite drug–disease and projected networks │ └── visualisations/ │ ├── drug_centrality_plot.py # Visualises drug centrality rankings │ ├── disease_projection_plot.py # Shows disease–disease co-occurrence via shared drugs │ └── candidate_insights.py # Highlights top repurposing candidates │ └── centrality_metrics.py # Calculates degree, betweenness, eigenvector centrality ├── data/ # JSON files with filtered trials per therapeutic area ├── graph/ # GraphML files for bipartite and projected networks ├── output/ │ ├── centrality_summary.csv # CSV file with centrality scores │ └── top_candidate_trials.csv # Disease-trial mappings for top drugs ├── plots/ # All generated figures └── README.md ``` </code> </pre>


## Methodology Overview

1. **Trial Filtering**  
   - Only **drug-based** interventional trials  
   - Status: `Completed`, `Recruiting`, `Active`  
   - Phases: `I`, `II`, `III`, `IV`

2. **Graph Construction**  
   - **Bipartite Network:** Drugs ↔ Diseases  
   - **Projections:**  
     - Drug–Drug: shared disease connections  
     - Disease–Disease: shared drugs

3. **Centrality Analysis**  
   - **Degree Centrality** – Exposure across diseases  
   - **Betweenness Centrality** – Bridge between disease clusters  
   - **Eigenvector Centrality** – Influence in core clusters

4. **Visualisation**  
   - Top drug hubs coloured and labelled  
   - Network edge labels denote shared entity counts  
   - Disease networks filtered by shared drug threshold

## Key Findings

- **Dexamethasone**, **Rituximab**, and **Prednisone** rank highly across all centrality metrics, suggesting strong repurposing potential.
- Disease clusters linked by common drugs reveal biological and therapeutic overlaps, e.g., inflammation, oncology, and autoimmunity.

## Plots Included

- Bipartite Drug–Disease Network
- Drug Centrality Highlights (Degree, Betweenness, Eigenvector)
- Disease–Disease Network (based on shared drugs)
- Centrality Summary Table (CSV)

## Requirements

- Python 3.8+
- `networkx`, `matplotlib`, `requests`, `pandas`

Install dependencies:

```bash
pip install -r requirements.txt



