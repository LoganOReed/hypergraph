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

    def reactions_to_genes(self, network):
        reactions_to_genes = {}
        for gene in network.genes:
            for reaction in gene.reaction.split(" "):
                reactions_to_genes[reaction] = gene.name.split(" ")

        return reactions_to_genes

    def get_well_ids_of_genes_in_network(self, genes_to_reactions_dict):
        # in this section, we are looking up the well_ids to tie back to the
        # drug_data in order to get the effect on the reactions
        genes_in_network = []
        for key in genes_to_reactions_dict.keys():
            genes_in_network.append(key)

        well_id_in_network_list = []
        for gene in genes_in_network:
            try:
                # gene is prefixed by mtu:, clone_ORF_lookup excludes this prefix
                index = clone_ORF_lookup.iloc[:, 1].tolist().index(gene[4:])
                well_id_in_network_list.append(clone_ORF_lookup.iloc[index, 0])
            except ValueError:
                pass

        return well_id_in_network_list


if __name__ == "__main__":
    absolutePath = os.path.abspath(__file__)
    fileDirectory = os.path.dirname(absolutePath)
    parentDirectory = os.path.dirname(fileDirectory)
    path_drug_data = os.path.join(
        parentDirectory, "input_files/Multidrug_6hr_Responses.xlsx"
    )

    path_KEGG = os.path.join(parentDirectory, "input_files/KEGG_data/")

    path_network = path_KEGG + "mtu01200.xml"

    path_cloneID_ORF = os.path.join(parentDirectory, "input_files/clone_to_orf.csv")

    excel = readXLSX(path_drug_data)
    drug_data = Drug_Data(excel)

    pathway_obj = read_KGML(path_network)
    network = KEGG_Network(pathway_obj)

    clone_ORF_lookup = read_cloneID_to_orf_table(path_cloneID_ORF)

    # print("######### reactions ##########")
    # for reaction in network.reactions:
    #     print(reaction)

    # print("########### genes ###############")
    # for gene in network.genes:
    #     print(gene)

    # print("######## misc ######")
    # gene = network.genes[10]
    # print(gene)
    # print(gene.reaction)

    d_u = Determining_Uberedges(drug_data)
    genes_to_reactions_dict = d_u.genes_to_reactions(network)
    reactions_to_genes_dict = d_u.reactions_to_genes(network)

    well_id_in_network_list = d_u.get_well_ids_of_genes_in_network(
        genes_to_reactions_dict)

    # this section will find the drug effect for each reaction in the list
    drug_combo = [10, 14, 17, 24]
    drug_table = drug_data.drug_table
    drug_effect_dict = {}
    for well_id in well_id_in_network_list:
        reverse_index = clone_ORF_lookup.iloc[:, 0].tolist().index(well_id)
        gene = clone_ORF_lookup.iloc[reverse_index, 1]

        index_of_well_in_drug_table = drug_table.iloc[:, 0].tolist().index(well_id)
        drug_effects = drug_table.iloc[index_of_well_in_drug_table, drug_combo]
        
        summed_drug_effect = sum(drug_effects)
        drug_effect_dict[gene] = summed_drug_effect

    # print("######### drug effect dictionary #########")
    # print(drug_effect_dict)

    # print("########### genes to reactions dictionary #########")
    # print(genes_to_reactions_dict)

    # this section will use the gene to reaction dict
    # and the drug effet dictionary to create a new dictionary
    # which contains reactions as keys and values of effect
    reaction_drug_effect_dict = {}
    for key in drug_effect_dict.keys():
        drug_effect = drug_effect_dict[key]
        new_key = "mtu:" + key
        reactions = genes_to_reactions_dict[new_key]
        for reaction in reactions:
            reaction_drug_effect_dict[reaction] = drug_effect

    print(reaction_drug_effect_dict)