import os
import networkx as nx
import pandas as pd

# Load the drug-drug network
graph_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "graph", "drug_network.graphml"))
G = nx.read_graphml(graph_path)

# Standardize node labels to lowercase
G = nx.relabel_nodes(G, {n: n.lower() for n in G.nodes()})

# Compute centrality measures
degree_centrality = nx.degree_centrality(G)
betweenness_centrality = nx.betweenness_centrality(G, weight="weight", normalized=True)
eigenvector_centrality = nx.eigenvector_centrality(G, weight="weight", max_iter=1000)

# Compile into a DataFrame
centrality_df = pd.DataFrame({
    "Drug": list(degree_centrality.keys()),
    "Degree Centrality": list(degree_centrality.values()),
    "Betweenness Centrality": list(betweenness_centrality.values()),
    "Eigenvector Centrality": list(eigenvector_centrality.values())
})

# Sort by Degree Centrality
centrality_df = centrality_df.sort_values(by="Degree Centrality", ascending=False).reset_index(drop=True)
filtered_df = centrality_df[~centrality_df["Drug"].str.contains("placebo", case=False)]
# Save or print summary
print(filtered_df.head(10))
