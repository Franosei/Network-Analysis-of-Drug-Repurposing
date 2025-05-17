import os
import networkx as nx
import matplotlib.pyplot as plt

# Paths
graph_dir = os.path.abspath(os.path.join(os.getcwd(), "..", "graph"))
drug_network_path = os.path.join(graph_dir, "drug_network.graphml")

# Load and relabel graph
G = nx.read_graphml(drug_network_path)
G = nx.relabel_nodes(G, {n: n.lower() for n in G.nodes()})

# Optional: Remove 'placebo' or similar non-drugs
G = G.subgraph([n for n in G.nodes if "placebo" not in n])

# Optional: Reduce size for performance (top 100 by degree)
top_k = 100
core_nodes = sorted(G.degree, key=lambda x: x[1], reverse=True)[:top_k]
G = G.subgraph([n for n, _ in core_nodes])

# Compute eigenvector centrality (ensure graph is connected or handle exception)
try:
    eigenvector = nx.eigenvector_centrality(G, max_iter=1000, weight="weight")
except nx.NetworkXException as e:
    print("Eigenvector centrality failed:", e)
    eigenvector = {n: 0 for n in G.nodes()}

# Identify top 20 drugs by eigenvector centrality
top_n = 20
top_drugs = sorted(eigenvector.items(), key=lambda x: x[1], reverse=True)[:top_n]
top_nodes = {n for n, _ in top_drugs}

# Node styling
node_sizes = [500 if n in top_nodes else 100 for n in G.nodes()]
node_colors = ["gold" if n in top_nodes else "skyblue" for n in G.nodes()]
labels = {n: n.title() for n in top_nodes}

# Layout
pos = nx.spring_layout(G, k=0.45, iterations=30, seed=42)

# Plot
plt.figure(figsize=(20, 15))
nx.draw_networkx_edges(G, pos, alpha=0.2, edge_color="gray")
nx.draw_networkx_nodes(G, pos, node_size=node_sizes, node_color=node_colors, alpha=0.9)
nx.draw_networkx_labels(G, pos, labels, font_size=9, font_color="black")

plt.title("Drug Network – Top 20 by Eigenvector Centrality", fontsize=16)
plt.axis("off")
plt.tight_layout()
plt.show()
