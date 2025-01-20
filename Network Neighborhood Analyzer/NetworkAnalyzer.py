"""
Author: P.A.R.R. Wijerathne
Date: 3/20/2024
Project 7: Network neighborhood analyzer
Inputs: A tsv file containing protein-protein interaction data,
        A tsv/text file containing protein information related to a known function.
Outputs: 
NetworkX graph object
The number of proteins and interactions in a given PPI network
The degree of a given protein in a given PPI network
The degree distribution as a histogram 
The number of proteins annotated to a given function in the immediate neighborhood
The distribution of annotated protein count in the immediate neighborhoods for all the proteins  
candidate genes for a single function using the Hishigaki algorithm
Hishigaki score distribution as a histogram.
Hishigaki score distribution as a boxplot

"""
import networkx as nx
import matplotlib.pyplot as plt
from collections import OrderedDict

class NetworkAnalyzer:
    def __init__(self, file_input):
        # Initialize with the file containing PPI network data
        self.file_input = file_input
        # Create the network graph
        self.G = self._create_network()

    def _create_network(self):
        # Internal method to create a NetworkX graph from the provided file
        G = nx.Graph()
        with open(self.file_input, "r") as tsv:
            all_set = set()
            for line in tsv:
                if line.startswith("#") or line.startswith(" "):
                    continue
                else:
                    record = line.upper().strip().split("\t")
                    node1, node2 = record[0], record[1]
                    G.add_edge(node1, node2, c_score=float(record[-1]))
                    all_set.update({node1, node2})
        return G

    def count_proteins_interactions(self):
        # Count the number of proteins and interactions in the network
        number_of_edges = len(self.G.edges)
        number_of_proteins = len(self.G.nodes)
        return f"Number of proteins: {number_of_proteins} and Number of Interactions: {number_of_edges} in the given PPI network."

    def calculate_degree(self, protein):
        # Calculate the degree of a given protein in the network
        return self.G.degree(protein)

    def plot_degree_distribution(self):
        # Plot the degree distribution of the network
        degree_list = [self.G.degree(node) for node in self.G.nodes]
        plt.figure(figsize=(8, 5))
        plt.hist(x=degree_list, bins=20, edgecolor='black')
        plt.xlabel("Number of interacting partners (Degree)")
        plt.ylabel("Number of proteins")
        plt.title("Degree Distribution of PPI Network")
        plt.grid(axis='y', linestyle='--', alpha=0.7)
        plt.show()

    def neighborhood_pcounter(self, protein_function_file):
        # Count the number of proteins annotated to a given function in the immediate neighborhood of each protein
        set_known = set()
        with open(protein_function_file, "r") as txt:
            for line in txt:
                recordx = line.upper().strip().split('\t')
                set_known.add(recordx[2])

        unknown_proteins = set(self.G.nodes) - set_known
        unknown_dict = {}
        for protein in unknown_proteins:
            neighbors = list(self.G.neighbors(protein))
            known_list = []
            for neighbor in neighbors:
                if neighbor in set_known:
                    known_list.append(neighbor)
            unknown_dict[protein] = len(known_list)

        return unknown_dict
    
    def neighborhood_pcounter_order(self,protein_function_file):
        # Order the proteins based on the count of proteins annotated to a given function in their immediate neighborhood
        dict1 = self.neighborhood_pcounter(protein_function_file)
        sortedDict = sorted(dict1.items(), key=lambda x:x[1], reverse=True)
        ordered = OrderedDict(sortedDict)
        # Extract top 5 proteins with highest counts
        best_hits = list(ordered.items())[0:5]
        return best_hits
    
    def annotated_protein_distribution(self, protein_function_file, remove_zeros):
        # Plot the distribution of annotated protein count in the immediate neighborhoods
        dict1 = self.neighborhood_pcounter(protein_function_file)
        protein_count = []
        dict_values = dict1.values()
        if remove_zeros:
            for value in dict_values:
                if value != 0:
                    protein_count.append(value)
        else:
            protein_count = dict_values

        plt.figure(figsize=(8, 5))
        plt.hist(x=protein_count, bins=20, edgecolor='black')
        plt.xlabel("Number of Annotated Protein Count in the Immediate Neighborhood")
        plt.ylabel("Number of Proteins")
        plt.title("Distribution of Annotated Protein Count ")
        plt.grid(axis='y', linestyle='--', alpha=0.7)
        plt.show()

    def protein_neighbors(self, protein_function_file,protein):
        # Analyze the immediate neighborhood of a given protein
        set_known = set()
        with open(protein_function_file, "r") as txt:
            for line in txt:
                recordx = line.upper().strip().split('\t')
                set_known.add(recordx[2])
        
        neighbors = list(self.G.neighbors(protein))
        for neighbor in neighbors:
            if neighbor in set_known:
                known_list = []
                for neighbor in neighbors:
                    if neighbor in set_known:
                        known_list.append(neighbor)
                          
        percentage = len(known_list)/len(neighbors)*100
        print(f"Number of annotated proteins in the immediate neighborhood of {protein} Protein: {len(known_list)}")
        print(f"{round(percentage,2)}% of proteins in the neighborhood are annotated")

    def hishigaki_method(self, protein_function_file):
        # Predict candidate genes for a single function using the Hishigaki algorithm
        hishigaki_scores = {}
        set_known = set()
        with open(protein_function_file, "r") as txt:
            for line in txt:
                recordx = line.upper().strip().split('\t')
                set_known.add(recordx[2])

        unknown_proteins = set(self.G.nodes) - set_known
        total_proteins = len(self.G.nodes)

        pc = 0
        for protein in self.G.nodes:
            if protein in set_known:
                pc = pc+1

        freq = (pc / total_proteins)

        for protein in unknown_proteins:
            neighbors = list(self.G.neighbors(protein))
            known_list = [neighbor for neighbor in neighbors if neighbor in set_known]
            nf = len(known_list)
            ef = freq * nf
            if ef == 0:
                hish_score = ((nf - ef) ** 2 / (ef + 0.00000001))
            else:
                hish_score = ((nf - ef) ** 2 / ef)
            hishigaki_scores[protein] = hish_score

        return(hishigaki_scores)
    
    def hishigaki_analyze(self, protein_function_file):
        # Analyze the Hishigaki scores
        hish_score_dict = self.hishigaki_method(protein_function_file)
        sortedDict = sorted(hish_score_dict.items(), key=lambda x: x[1], reverse=True)
        ordered = OrderedDict(sortedDict)

        best_hit = list(ordered.items())[0:5]
        return best_hit

    def plot_hishigaki_score_distribution(self, protein_function_file):
        # Plot the distribution of Hishigaki scores
        hish_score_dict = self.hishigaki_method(protein_function_file)
        hishigaki_scores = hish_score_dict.values()

        plt.hist(hishigaki_scores, bins=20, edgecolor='black')
        plt.xlabel("Hishigaki Score")
        plt.ylabel("Frequency")
        plt.title("Hishigaki Score Distribution")
        plt.grid(axis='y', linestyle='--', alpha=0.7)
        plt.show()

    def hishigaki_box_plot(self, protein_function_file):
        # Plot the distribution of Hishigaki scores(Box Plot)
        hish_score_dict = self.hishigaki_method(protein_function_file)
        hishigaki_scores = list(hish_score_dict.values())

        plt.boxplot(hishigaki_scores)
        plt.xlabel("Hishigaki Score")
        plt.ylabel("Score")
        plt.title("Distribution of Hishigaki Scores")
        plt.grid(axis='y', linestyle='--', alpha=0.7)
        plt.show()

