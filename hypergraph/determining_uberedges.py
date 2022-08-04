import pandas as pd
import os
from drug_data import Drug_Data
from kegg_network import KEGG_Network
from fileIO import readXLSX, read_KGML, read_cloneID_to_orf_table


class Determining_Uberedges:
    def __init__(self, drug_data):
        self.drug_data = drug_data

    def genes_to_reactions(self, network):
        # creates a dictionary where key is gene name
        # and values are lists of reactions that are effected
        # Note: some genes in network.genes have more than one name, we split
        # based on ' ' to have keys that correspond to exactly one name
        genes_to_reactions = {}
        for gene in network.genes:
            for name in gene.name.split(" "):
                genes_to_reactions[name] = gene.reaction.split(" ")

        return genes_to_reactions


if __name__ == "__main__":
    absolutePath = os.path.abspath(__file__)
    fileDirectory = os.path.dirname(absolutePath)
    parentDirectory = os.path.dirname(fileDirectory)
    path_drug_data = os.path.join(
        parentDirectory, "input_files/Multidrug_6hr_Responses_trimmed.xlsx"
    )

    path_KEGG = os.path.join(parentDirectory, "input_files/KEGG_data/")

    path_network = path_KEGG + "mtu01200.xml"

    path_cloneID_ORF = os.path.join(parentDirectory, "input_files/clone_to_orf.csv")

    excel = readXLSX(path_drug_data)
    drug_data = Drug_Data(excel)

    pathway_obj = read_KGML(path_network)
    network = KEGG_Network(pathway_obj)

    clone_ORF_lookup = read_cloneID_to_orf_table(path_cloneID_ORF)

    print("######### reactions ##########")
    for reaction in network.reactions:
        print(reaction)

    print("########### genes ###############")
    for gene in network.genes:
        print(gene)

    print("######## misc ######")
    gene = network.genes[10]
    print(gene)
    print(gene.reaction)

    d_u = Determining_Uberedges(drug_data)
    test = d_u.genes_to_reactions(network)

    print(test)