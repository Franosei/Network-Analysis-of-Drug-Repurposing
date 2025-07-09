import os
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt


def clustered_bipartite_layout(G, left_nodes, right_nodes, x_gap=6, y_gap=1.5):
    pos = {}
    max_len = max(len(left_nodes), len(right_nodes))
    left_y = np.linspace(0, y_gap * (max_len - 1), len(left_nodes))
    right_y = np.linspace(0, y_gap * (max_len - 1), len(right_nodes))
    for i, node in enumerate(sorted(left_nodes)):
        pos[node] = (0, left_y[i])
    for i, node in enumerate(sorted(right_nodes)):
        pos[node] = (x_gap, right_y[i])
    return pos

def plot_top50_clustered_colored(
    G_full,
    output_path="figures/bipartite_top50_clustered_colored.png",
    top_n=50,
    figsize=(20, 18),
    dpi=600
):
    # Identify node types
    drug_nodes = [n for n, d in G_full.nodes(data=True) if d.get("type") == "drug"]
    disease_nodes = [n for n, d in G_full.nodes(data=True) if d.get("type") == "disease"]

    # Select top-N by degree
    degree_dict = dict(G_full.degree())
    top_drugs = sorted(drug_nodes, key=lambda x: degree_dict[x], reverse=True)[:top_n]
    top_diseases = sorted(disease_nodes, key=lambda x: degree_dict[x], reverse=True)[:top_n]

    # Create subgraph
    sub_nodes = top_drugs + top_diseases
    G = G_full.subgraph(sub_nodes).copy()

    # Layout
    pos = clustered_bipartite_layout(G, top_drugs, top_diseases)

    # Assign unique color to each drug
    cmap = plt.colormaps.get_cmap('tab20').resampled(top_n)
    drug_color_map = {drug: cmap(i) for i, drug in enumerate(top_drugs)}

    # Node coloring
    node_colors = [drug_color_map[n] if n in drug_color_map else '#e31a1c' for n in G.nodes()]
    node_sizes = [600 for _ in G.nodes()]

    # Plot
    plt.figure(figsize=figsize)

    # Draw colored edges
    for drug in top_drugs:
        targets = list(G.neighbors(drug))
        edges = [(drug, t) for t in targets]
        nx.draw_networkx_edges(
            G,
            pos,
            edgelist=edges,
            edge_color=[drug_color_map[drug]] * len(edges),
            alpha=0.4,
            width=0.9
        )

    # Draw nodes and labels
    nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=node_sizes, alpha=0.95)
    nx.draw_networkx_labels(G, pos, font_size=10)

    # Legend (top 5 drugs + disease)
    legend_drugs = list(top_drugs[:5])
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', label=drug, markerfacecolor=drug_color_map[drug], markersize=8)
        for drug in legend_drugs
    ]
    legend_elements.append(Line2D([0], [0], marker='o', color='w', label='Disease', markerfacecolor='#e31a1c', markersize=8))
    plt.legend(handles=legend_elements, loc='upper right', fontsize=8, title="Top Drugs")

    plt.title("Top 50 Drug–Disease Bipartite Graph (Clustered Layout, Drug-Colored Edges)", fontsize=15)
    plt.axis('off')
    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi)
    print(f"Top-50 clustered plot saved to: {os.path.abspath(output_path)}")

def plot_full_auto_layout(
    G,
    output_path="figures/bipartite_full_auto.png",
    layout="spring",
    figsize=(18, 14),
    dpi=600
):
    drug_nodes = [n for n, d in G.nodes(data=True) if d.get("type") == "drug"]
    disease_nodes = [n for n, d in G.nodes(data=True) if d.get("type") == "disease"]

    pos = {
        "spring": nx.spring_layout,
        "kamada_kawai": nx.kamada_kawai_layout,
        "spectral": nx.spectral_layout
    }.get(layout, nx.spring_layout)(G, seed=42)

    node_colors = ['#1f78b4' if n in drug_nodes else '#e31a1c' for n in G.nodes()]
    node_sizes = [250 for _ in G.nodes()]

    plt.figure(figsize=figsize)
    nx.draw_networkx_edges(G, pos, alpha=0.05, width=0.4)
    nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=node_sizes, alpha=0.85)
    nx.draw_networkx_labels(G, pos, font_size=8)

    legend = [
        Line2D([0], [0], marker='o', color='w', label='Drug', markerfacecolor='#1f78b4', markersize=7),
        Line2D([0], [0], marker='o', color='w', label='Disease', markerfacecolor='#e31a1c', markersize=7)
    ]
    plt.legend(handles=legend, loc='upper right', fontsize=8)
    plt.title("Full Drug–Disease Bipartite Network (Spring Layout)", fontsize=14)
    plt.axis('off')
    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi)
    print(f"Full auto-layout plot saved to: {os.path.abspath(output_path)}")

if __name__ == "__main__":
    os.makedirs("figures", exist_ok=True)
    graph_path = "graph/bipartite.graphml"
    G = nx.read_graphml(graph_path)

    # Plot 1: Clustered top 50 with colored drug edges
    plot_top50_clustered_colored(G)

    # Plot 2: Full graph with auto layout (unchanged)
    plot_full_auto_layout(G)
