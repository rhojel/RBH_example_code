"""
This file creates a graph visualization of a network of proteins. The proteins are connected based on 
a significance of the connection (calculated in networkgen.py). The proteins are colored according 
to their indications.

Usage: plot_graph.py <network_edgelist.txt> <protein_nodes.csv> <file_path_to_output_figure.png>

"""

from chemoUtils import readNodeFile
import sys
import networkx as nx
import matplotlib.pyplot as plt

def main():
    """
    Main function that reads in the input files into a graph, generates a color mapping, and graphs
    the protein network. The graph is then saved to the directory included in the user input.

    """
    node_dict, color_dict = readNodeFile(sys.argv[2],True)

    H = nx.read_edgelist(sys.argv[1])
    G = nx.relabel_nodes(H, node_dict)
    G.add_nodes_from(node_dict[k] for k in node_dict)   # adds any nodes that aren't in the edge list 

    color_map = []
    for node in G:
        if color_dict[node] == 'bp':
            color_map.append('red')
        elif color_dict[node] == 'bp;cholesterol':
            color_map.append('green')
        elif color_dict[node] == 'bp;cholesterol;diabetes':
            color_map.append('blue')
        else:
            color_map.append('purple')

    plt.figure(figsize=(8, 8), dpi=150)
    nx.draw_networkx(G,node_color=color_map, with_labels=True, pos = nx.spring_layout(G,k=0.8,iterations=50), alpha = 0.8)
    plt.axis('off')
    plt.savefig(sys.argv[3])

if __name__ == '__main__':
    main()
