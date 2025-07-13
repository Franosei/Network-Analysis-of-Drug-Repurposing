import os
import json
import networkx as nx
import pandas as pd
import numpy as np
from scipy.sparse import csr_matrix, identity
from scipy.sparse.linalg import inv


class InterpretableGraphFeatureBuilder:
    """
    Build a bipartite drug–disease graph and compute five interpretable, information-rich graph features
    for both known and unknown drug–disease pairs:

    1. Graph Distance to Known Indications (inverse shortest path)
    2. Random Walk with Restart (personalized PageRank from drug to disease)
    3. Structural Likelihood (centrality(drug) × centrality(disease))
    4. Preferential Attachment (degree(drug) × degree(disease))
    5. Katz Similarity (true matrix-based similarity from all damped paths)

    Outputs:
    - graph_features_known.csv: All known trial pairs (label = 1)
    - graph_features_unknown.csv: All novel/unseen drug–disease pairs (label = 0)
    - bipartite.graphml: GraphML file of the bipartite network for visualization
    """

    def __init__(self, input_file="processed_data/condition_drug_pairs.json",
                 known_output="graph/graph_features_known.csv",
                 unknown_output="graph/graph_features_unknown.csv"):
        self.input_file = input_file
        self.known_output = known_output
        self.unknown_output = unknown_output
        os.makedirs(os.path.dirname(known_output), exist_ok=True)

        self.graph = nx.Graph()
        self.drugs = set()
        self.diseases = set()
        self.known_indications = {}

    def build_bipartite_graph(self):
        with open(self.input_file, "r", encoding="utf-8") as f:
            data = json.load(f)

        for entry in data:
            drug = entry["intervention"].strip().lower()
            disease = entry["condition"].strip().lower()

            if "placebo" in drug:
                continue

            self.graph.add_node(drug, bipartite="drug")
            self.graph.add_node(disease, bipartite="disease")
            self.graph.add_edge(drug, disease)

            self.drugs.add(drug)
            self.diseases.add(disease)

            if drug not in self.known_indications:
                self.known_indications[drug] = set()
            self.known_indications[drug].add(disease)

        print(f"Graph: {len(self.drugs)} drugs, {len(self.diseases)} diseases, {len(self.graph.edges())} edges")

        # Save the graph to GraphML
        nx.write_graphml(self.graph, "graph/bipartite.graphml")
        print("Saved graph to graph/bipartite.graphml")

    def compute_katz_similarity_matrix(self, alpha=0.005):
        nodes = list(self.graph.nodes())
        node_index = {n: i for i, n in enumerate(nodes)}
        index_node = {i: n for n, i in node_index.items()}

        A = nx.to_scipy_sparse_array(self.graph, nodelist=nodes, format='csr')
        I = identity(A.shape[0], format='csr')
        try:
            K = inv(I - alpha * A) - I
            K = K.toarray()
        except Exception as e:
            print(f"⚠️ Katz matrix computation failed: {e}")
            return {}

        katz_scores = {}
        for i, drug in enumerate(nodes):
            if drug not in self.drugs:
                continue
            katz_scores[drug] = {}
            for j, disease in enumerate(nodes):
                if disease not in self.diseases:
                    continue
                katz_scores[drug][disease] = K[i, j]
        return katz_scores

    def compute_features_for_pair(self, drug, disease, degree, eigen, between,
                                  degree_disease, eigen_disease, between_disease,
                                  pageranks, katz_matrix):
        # Feature 1: Graph Distance to Known Indications
        try:
            min_dist = min(
                nx.shortest_path_length(self.graph, source=disease, target=ind)
                for ind in self.known_indications.get(drug, [])
                if nx.has_path(self.graph, disease, ind)
            )
            graph_distance_score = 1 / (1 + min_dist)
        except ValueError:
            graph_distance_score = 0.0

        # Feature 2: Random Walk with Restart (RWR)
        rwr_score = pageranks[drug].get(disease, 0.0)

        # Feature 3: Structural Likelihood
        centrality_drug = (
            0.33 * degree.get(drug, 0) +
            0.33 * eigen.get(drug, 0) +
            0.34 * between.get(drug, 0)
        )
        centrality_disease = (
            0.33 * degree_disease.get(disease, 0) +
            0.33 * eigen_disease.get(disease, 0) +
            0.34 * between_disease.get(disease, 0)
        )
        structural_likelihood = round((1 + centrality_drug) * (1 + centrality_disease), 4)

        # Feature 4: Preferential Attachment
        pa_score = degree.get(drug, 0) * degree_disease.get(disease, 0)

        # Feature 5: Katz Similarity (Matrix-based)
        katz_score = katz_matrix.get(drug, {}).get(disease, 0.0)

        return {
            "Drug": drug,
            "Disease": disease,
            "GraphDistanceToIndication": round(graph_distance_score, 4),
            "RandomWalkScore": round(rwr_score, 6),
            "StructuralLikelihood": structural_likelihood,
            "PreferentialAttachment": round(pa_score, 6),
            "KatzSimilarity": round(katz_score, 6)
        }

    def compute_all_features(self):
        known_rows = []
        unknown_rows = []

        # Centralities
        degree = nx.degree_centrality(self.graph)
        eigen = nx.eigenvector_centrality(self.graph, max_iter=1000)
        between = nx.betweenness_centrality(self.graph)
        degree_disease = {n: degree[n] for n in self.diseases}
        eigen_disease = {n: eigen[n] for n in self.diseases}
        between_disease = {n: between[n] for n in self.diseases}

        # RWR cache
        pageranks = {
            drug: nx.pagerank(self.graph, alpha=0.85, personalization={drug: 1.0})
            for drug in self.drugs
        }

        # Katz Similarity Matrix
        print("Computing matrix-based Katz similarity...")
        katz_matrix = self.compute_katz_similarity_matrix(alpha=0.005)
        print("Katz matrix computed.")

        # Loop over all drug–disease pairs
        for drug in self.drugs:
            for disease in self.diseases:
                row = self.compute_features_for_pair(
                    drug, disease,
                    degree, eigen, between,
                    degree_disease, eigen_disease, between_disease,
                    pageranks, katz_matrix
                )
                if disease in self.known_indications.get(drug, set()):
                    row["Label"] = 1
                    known_rows.append(row)
                else:
                    row["Label"] = 0
                    unknown_rows.append(row)

        # Save to CSV
        pd.DataFrame(known_rows).to_csv(self.known_output, index=False)
        pd.DataFrame(unknown_rows).to_csv(self.unknown_output, index=False)
        print(f"Saved known features to: {self.known_output}")
        print(f"Saved unknown features to: {self.unknown_output}")
        print(pd.DataFrame(unknown_rows).head(5))


if __name__ == "__main__":
    builder = InterpretableGraphFeatureBuilder()
    builder.build_bipartite_graph()
    builder.compute_all_features()
