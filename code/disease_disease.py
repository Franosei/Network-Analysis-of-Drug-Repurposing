import os
import networkx as nx
import matplotlib.pyplot as plt

# Load the bipartite graph
graph_dir = os.path.abspath(os.path.join(os.getcwd(), "..", "graph"))
bipartite_path = os.path.join(graph_dir, "bipartite.graphml")

B = nx.read_graphml(bipartite_path)
B = nx.relabel_nodes(B, {n: n.lower() for n in B.nodes()})

# Extract disease nodes
disease_nodes = {n for n, d in B.nodes(data=True) if d.get("bipartite") == "disease"}

# Create weighted disease–disease projection
disease_graph = nx.bipartite.weighted_projected_graph(B, disease_nodes)

# Filter edges with weight > 2 (i.e., shared drugs > 2)
filtered_edges = [(u, v, d) for u, v, d in disease_graph.edges(data=True) if d['weight'] > 5]
filtered_graph = nx.Graph()
filtered_graph.add_edges_from([(u, v, {'weight': d['weight']}) for u, v, d in filtered_edges])

# Layout and styles
pos = nx.spring_layout(filtered_graph, seed=42, k=1.2)
node_size = 2200
node_color = "gold"
label_color = "black"

# Draw nodes and edges
plt.figure(figsize=(10, 8))
nx.draw_networkx_nodes(filtered_graph, pos, node_color=node_color, node_size=node_size)
nx.draw_networkx_edges(filtered_graph, pos, edge_color="gray", width=1.5)

# Add labels
nx.draw_networkx_labels(filtered_graph, pos, font_size=10, font_color=label_color)

# Edge labels = number of shared drugs
edge_labels = nx.get_edge_attributes(filtered_graph, 'weight')
nx.draw_networkx_edge_labels(filtered_graph, pos, edge_labels=edge_labels, font_size=9)

plt.title("Disease–Disease Network (Shared Drugs > 5)", fontsize=16)
plt.axis("off")
plt.tight_layout()
plt.show()
