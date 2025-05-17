import os
import csv
import networkx as nx

def compute_and_save_drug_centrality(graph_path, output_csv):
    # Load graph
    G = nx.read_graphml(graph_path)
    G = nx.relabel_nodes(G, {n: n.lower() for n in G.nodes()})

    # Identify drug nodes (excluding placebo)
    drug_nodes = {n for n, d in G.nodes(data=True) if d.get("bipartite") == "drug" and "placebo" not in n}

    # Compute degree centrality
    degree_centrality = {node: G.degree(node) for node in drug_nodes}
    sorted_centrality = sorted(degree_centrality.items(), key=lambda x: x[1], reverse=True)

    # Save to CSV
    with open(output_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["Drug", "DegreeCentrality"])
        writer.writerows(sorted_centrality)

    print(f"Centrality rankings saved to: {output_csv}")

if __name__ == "__main__":
    base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "graph"))
    graph_file = os.path.join(base_dir, "bipartite.graphml")
    output_file = os.path.join(base_dir, "drug_centrality.csv")

    compute_and_save_drug_centrality(graph_file, output_file)
