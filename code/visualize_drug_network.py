import os
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np


class DrugDiseaseGraphVisualizer:
    """
    Visualize a drug–disease bipartite graph from a saved GraphML file.
    It provides:
    - Full graph visualization with layout options
    - Top-N clustered bipartite layout visualization with colored drug edges
    """

    def __init__(self, graphml_path="graph/bipartite.graphml", output_dir="figures"):
        self.graphml_path = graphml_path
        self.output_dir = output_dir
        os.makedirs(self.output_dir, exist_ok=True)
        self.G = nx.read_graphml(graphml_path)
        self._ensure_node_types()

    def _ensure_node_types(self):
        """Ensure all nodes are tagged with type='drug' or 'disease' based on bipartite attribute"""
        for node, attr in self.G.nodes(data=True):
            bipartite_value = attr.get("bipartite", "")
            self.G.nodes[node]["type"] = bipartite_value

    def _clustered_bipartite_layout(self, left_nodes, right_nodes, x_gap=6, y_gap=1.5):
        """Manual layout for bipartite clusters: drugs left, diseases right"""
        pos = {}
        max_len = max(len(left_nodes), len(right_nodes))
        left_y = np.linspace(0, y_gap * (max_len - 1), len(left_nodes))
        right_y = np.linspace(0, y_gap * (max_len - 1), len(right_nodes))
        for i, node in enumerate(sorted(left_nodes)):
            pos[node] = (0, left_y[i])
        for i, node in enumerate(sorted(right_nodes)):
            pos[node] = (x_gap, right_y[i])
        return pos

    def plot_topN_clustered_colored(self, top_n=50, filename="bipartite_topN_clustered_colored.png"):
        """Plot a top-N bipartite graph with clustered layout and drug-colored edges"""
        drug_nodes = [n for n, d in self.G.nodes(data=True) if d["type"] == "drug"]
        disease_nodes = [n for n, d in self.G.nodes(data=True) if d["type"] == "disease"]

        degree_dict = dict(self.G.degree())
        top_drugs = sorted(drug_nodes, key=lambda x: degree_dict.get(x, 0), reverse=True)[:top_n]
        top_diseases = sorted(disease_nodes, key=lambda x: degree_dict.get(x, 0), reverse=True)[:top_n]

        sub_nodes = top_drugs + top_diseases
        G_sub = self.G.subgraph(sub_nodes).copy()
        pos = self._clustered_bipartite_layout(top_drugs, top_diseases)

        cmap = plt.colormaps.get_cmap('tab20').resampled(top_n)
        drug_color_map = {drug: cmap(i) for i, drug in enumerate(top_drugs)}

        node_colors = [drug_color_map[n] if n in drug_color_map else '#e31a1c' for n in G_sub.nodes()]
        node_sizes = [600 for _ in G_sub.nodes()]

        plt.figure(figsize=(20, 18))

        for drug in top_drugs:
            targets = list(G_sub.neighbors(drug))
            edges = [(drug, t) for t in targets]
            nx.draw_networkx_edges(
                G_sub, pos, edgelist=edges,
                edge_color=[drug_color_map[drug]] * len(edges),
                alpha=0.4, width=0.9
            )

        nx.draw_networkx_nodes(G_sub, pos, node_color=node_colors, node_size=node_sizes, alpha=0.95)
        nx.draw_networkx_labels(G_sub, pos, font_size=10)

        legend_elements = [
            Line2D([0], [0], marker='o', color='w', label=drug,
                   markerfacecolor=drug_color_map[drug], markersize=8)
            for drug in top_drugs[:5]
        ]
        legend_elements.append(Line2D([0], [0], marker='o', color='w', label='Disease',
                                      markerfacecolor='#e31a1c', markersize=8))
        plt.legend(handles=legend_elements, loc='upper right', fontsize=8, title="Top Drugs")

        plt.title(f"Top {top_n} Drug–Disease Bipartite Graph (Clustered Layout, Drug-Colored Edges)", fontsize=15)
        plt.axis('off')
        plt.tight_layout()
        save_path = os.path.join(self.output_dir, filename)
        plt.savefig(save_path, dpi=600)
        plt.close()
        print(f"Saved clustered top-{top_n} plot to: {save_path}")

    def plot_full_auto_layout(self, layout="spring", filename="bipartite_full_auto.png"):
        """Plot the full bipartite graph with auto layout (spring, kamada_kawai, spectral)"""
        layout_fn = {
            "spring": nx.spring_layout,
            "kamada_kawai": nx.kamada_kawai_layout,
            "spectral": nx.spectral_layout
        }.get(layout, nx.spring_layout)

        pos = layout_fn(self.G, seed=42)

        drug_nodes = [n for n, d in self.G.nodes(data=True) if d["type"] == "drug"]
        disease_nodes = [n for n, d in self.G.nodes(data=True) if d["type"] == "disease"]

        node_colors = ['#1f78b4' if n in drug_nodes else '#e31a1c' for n in self.G.nodes()]
        node_sizes = [250 for _ in self.G.nodes()]

        plt.figure(figsize=(18, 14))
        nx.draw_networkx_edges(self.G, pos, alpha=0.05, width=0.4)
        nx.draw_networkx_nodes(self.G, pos, node_color=node_colors, node_size=node_sizes, alpha=0.85)
        nx.draw_networkx_labels(self.G, pos, font_size=8)

        legend = [
            Line2D([0], [0], marker='o', color='w', label='Drug', markerfacecolor='#1f78b4', markersize=7),
            Line2D([0], [0], marker='o', color='w', label='Disease', markerfacecolor='#e31a1c', markersize=7)
        ]
        plt.legend(handles=legend, loc='upper right', fontsize=8)
        plt.title("Full Drug–Disease Bipartite Network (Auto Layout)", fontsize=14)
        plt.axis('off')
        plt.tight_layout()
        save_path = os.path.join(self.output_dir, filename)
        plt.savefig(save_path, dpi=600)
        plt.close()
        print(f"Saved full auto-layout plot to: {save_path}")


if __name__ == "__main__":
    visualizer = DrugDiseaseGraphVisualizer(graphml_path="graph/bipartite.graphml")
    visualizer.plot_topN_clustered_colored(top_n=50)
    visualizer.plot_full_auto_layout(layout="spring")
