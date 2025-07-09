import os
import json
import networkx as nx
import pandas as pd

class DrugRepurposingNetworkBuilder:
    def __init__(self, input_file="processed_data/condition_drug_pairs.json", graph_dir="graph"):
        self.input_file = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", input_file))
        self.graph_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", graph_dir))
        os.makedirs(self.graph_dir, exist_ok=True)

        self.bipartite_graph = nx.Graph()
        self.drug_nodes = set()
        self.disease_nodes = set()
        self.placebo_skipped = 0

    def build_bipartite_network(self):
        with open(self.input_file, "r", encoding="utf-8") as f:
            pairs = json.load(f)

        for item in pairs:
            drug_raw = item["intervention"].strip()
            condition_raw = item["condition"].strip()

            if "placebo" in drug_raw.lower():
                self.placebo_skipped += 1
                continue

            drug = drug_raw.lower()
            condition = condition_raw.lower()

            self.bipartite_graph.add_node(drug, bipartite="drug", type="drug")
            self.bipartite_graph.add_node(condition, bipartite="disease", type="disease")
            self.bipartite_graph.add_edge(drug, condition)

            self.drug_nodes.add(drug)
            self.disease_nodes.add(condition)

        print(f"\n Bipartite Graph Created:")
        print(f"  Drugs:    {len(self.drug_nodes)}")
        print(f"  Diseases: {len(self.disease_nodes)}")
        print(f"  Edges:    {self.bipartite_graph.number_of_edges()}")
        print(f"  Placebo interventions skipped: {self.placebo_skipped}")

    def project_drug_network(self):
        projected = nx.bipartite.weighted_projected_graph(self.bipartite_graph, self.drug_nodes)
        for node in projected.nodes():
            projected.nodes[node]["type"] = "drug"

        print(f"\n Projected Drug Network:")
        print(f"  Drugs: {len(projected.nodes)}")
        print(f"  Edges: {len(projected.edges)}")
        return projected

    def compute_and_save_centralities(self, G, filename="drug_centrality.csv"):
        print("\nComputing centrality metrics...")

        degree_centrality = nx.degree_centrality(G)
        eigen_centrality = nx.eigenvector_centrality(G, max_iter=1000)
        betweenness_centrality = nx.betweenness_centrality(G)
        clustering = nx.clustering(G)
        pagerank = nx.pagerank(G, alpha=0.85)

        rows = []
        for node in G.nodes():
            rows.append({
                "Drug": node,
                "DegreeCentrality": degree_centrality.get(node, 0),
                "EigenvectorCentrality": eigen_centrality.get(node, 0),
                "BetweennessCentrality": betweenness_centrality.get(node, 0),
                "ClusteringCoefficient": clustering.get(node, 0),
                "RandomWalkCentrality": pagerank.get(node, 0)
            })

        df = pd.DataFrame(rows)
        output_path = os.path.join(self.graph_dir, filename)
        df.to_csv(output_path, index=False)
        print(f"Saved drug centrality CSV to: {output_path}")

    def save_networks(self, bipartite_filename="bipartite.graphml", drug_filename="drug_network.graphml"):
        bipartite_path = os.path.join(self.graph_dir, bipartite_filename)
        drug_path = os.path.join(self.graph_dir, drug_filename)

        nx.write_graphml(self.bipartite_graph, bipartite_path)
        print(f"\nSaved bipartite graph to: {bipartite_path}")

        drug_graph = self.project_drug_network()
        nx.write_graphml(drug_graph, drug_path)
        print(f"Saved drug–drug graph to: {drug_path}")

        self.compute_and_save_centralities(drug_graph)

if __name__ == "__main__":
    builder = DrugRepurposingNetworkBuilder()
    builder.build_bipartite_network()
    builder.save_networks()
