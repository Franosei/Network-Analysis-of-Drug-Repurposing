import os
import json
import networkx as nx

class DrugRepurposingNetworkBuilder:
    def __init__(self, data_dir="data", graph_dir="graph"):
        # Paths for input data and output graphs
        self.data_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", data_dir))
        self.graph_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", graph_dir))
        os.makedirs(self.graph_dir, exist_ok=True)

        self.bipartite_graph = nx.Graph()
        self.drug_nodes = set()
        self.disease_nodes = set()

    def load_data(self):
        """
        Load all JSON trial files from the data directory.
        """
        trial_files = [f for f in os.listdir(self.data_dir) if f.endswith(".json")]
        print(f"Found {len(trial_files)} therapeutic area files.")
        return trial_files

    def build_bipartite_network(self):
        """
        Build a bipartite graph linking drugs and diseases.
        Nodes: drugs and diseases
        Edges: a drug tested for a disease
        """
        for file in self.load_data():
            disease_name = os.path.splitext(file)[0].replace("_", " ")
            path = os.path.join(self.data_dir, file)

            with open(path, "r", encoding="utf-8") as f:
                trials = json.load(f)

            for trial in trials:
                interventions = trial.get("interventions", [])
                types = trial.get("interventionTypes", [])

                for drug, dtype in zip(interventions, types):
                    if dtype.upper() != "DRUG":
                        continue

                    drug = drug.strip().lower()
                    disease = disease_name.strip().lower()

                    self.bipartite_graph.add_node(drug, bipartite="drug")
                    self.bipartite_graph.add_node(disease, bipartite="disease")
                    self.bipartite_graph.add_edge(drug, disease)

                    self.drug_nodes.add(drug)
                    self.disease_nodes.add(disease)

        print(f"Bipartite Graph: {len(self.drug_nodes)} drugs, {len(self.disease_nodes)} diseases, {self.bipartite_graph.number_of_edges()} connections")

    def project_drug_network(self):
        """
        Project the bipartite graph into a drug–drug network:
        - Drugs are connected if they share a disease.
        - Edge weight = number of shared diseases.
        """
        projected = nx.bipartite.weighted_projected_graph(self.bipartite_graph, self.drug_nodes)
        print(f"Projected Drug Network: {len(projected.nodes)} drugs, {len(projected.edges)} edges")
        return projected

    def save_networks(self, bipartite_filename="bipartite.graphml", drug_filename="drug_network.graphml"):
        """
        Save both networks to the 'graph' folder in .graphml format.
        """
        bipartite_path = os.path.join(self.graph_dir, bipartite_filename)
        drug_path = os.path.join(self.graph_dir, drug_filename)

        nx.write_graphml(self.bipartite_graph, bipartite_path)
        print(f"Bipartite graph saved to {bipartite_path}")

        drug_graph = self.project_drug_network()
        nx.write_graphml(drug_graph, drug_path)
        print(f"Drug-drug graph saved to {drug_path}")

if __name__ == "__main__":
    builder = DrugRepurposingNetworkBuilder()
    builder.build_bipartite_network()
    builder.save_networks()
