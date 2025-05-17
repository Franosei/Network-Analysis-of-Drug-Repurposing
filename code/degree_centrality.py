import os
import networkx as nx
import matplotlib.pyplot as plt

# Load the bipartite graph
graph_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "graph", "bipartite.graphml"))
G = nx.read_graphml(graph_path)

# Convert all node names to lowercase for uniformity
G = nx.relabel_nodes(G, {n: n.lower() for n in G.nodes()})

# Filter out placebo and get drug nodes only
drug_nodes = {n for n, d in G.nodes(data=True) if d.get("bipartite") == "drug" and "placebo" not in n}

# Get top 10 drugs by degree
top_drugs = sorted(drug_nodes, key=lambda n: G.degree(n), reverse=True)[:20]

# Prepare node colors and sizes
node_color = []
node_size = []
node_labels = {}

for node in G.nodes():
    node_data = G.nodes[node]
    if node in top_drugs:
        node_color.append("gold")
        node_size.append(800)
        node_labels[node] = node
    elif node_data.get("bipartite") == "drug":
        node_color.append("green")
        node_size.append(20)
    else:
        node_color.append("red")  # disease nodes
        node_size.append(50)

# Position layout
pos = nx.spring_layout(G, k=0.40, iterations=50, seed=42)

# Plotting
plt.figure(figsize=(20, 15))
nx.draw_networkx_edges(G, pos, alpha=0.4, edge_color="gray")
nx.draw_networkx_nodes(G, pos, node_color=node_color, node_size=node_size)
nx.draw_networkx_labels(G, pos, labels=node_labels, font_size=8, font_color="black")

plt.title("Top 20 Drug Candidates and Their Associated Diseases", fontsize=14)
plt.axis("off")
plt.tight_layout()
plt.show()
