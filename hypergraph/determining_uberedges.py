import pandas as pd
import os
from drug_data import Drug_Data
from kegg_network import KEGG_Network
from fileIO import readXLSX, read_KGML, read_cloneID_to_orf_table


class Determining_Uberedges:

    """Class for determining which edges should be uber edges
    
    Attributes
    ----------
    clone_ORF_lookup : dataframe
        dataframe that maps well ids to TB genes
    drug_data : Drug Data
        Custom Drug Data object that holds DE data
    gene_drug_effect_dict : dictionary
        dictionary to map gene to drug effect
    genes_to_reactions_dict : dictionary
        dictionary to map genes to corresponding reactions
    network : KEGG Network
        Custom KEGG Network object, holds general network info
    reaction_drug_effect_dict : dictionary
        dictionary to map reaction name to drug effect
    reactions_to_genes_dict : dictionary
        dictionary to map reactions to genes
    well_id_in_network_list : list
        gives the index of genes in network iin the well id
    """

    def __init__(self, drug_data, clone_ORF_lookup, network):
        self.drug_data = drug_data
        self.network = network
        self.clone_ORF_lookup = clone_ORF_lookup

        self.genes_to_reactions_dict = self.genes_to_reactions()
        self.reactions_to_genes_dict = self.reactions_to_genes()

        self.well_id_in_network_list = self.generate_well_ids_of_genes_in_network()

    def genes_to_reactions(self):
        # creates a dictionary where key is gene name
        # and values are lists of reactions that are effected
        # Note: some genes in network.genes have more than one name, we split
        # based on ' ' to have keys that correspond to exactly one name
        genes_to_reactions = {}
        for gene in self.network.genes:
            for name in gene.name.split(" "):
                genes_to_reactions[name] = gene.reaction.split(" ")

        return genes_to_reactions

    def reactions_to_genes(self):
        reactions_to_genes = {}
        for gene in self.network.genes:
            for reaction in gene.reaction.split(" "):
                reactions_to_genes[reaction] = gene.name.split(" ")

        return reactions_to_genes

    def generate_well_ids_of_genes_in_network(self):
        # in this section, we are looking up the well_ids to tie back to the
        # drug_data in order to get the effect on the reactions
        genes_in_network = []
        for key in self.genes_to_reactions_dict.keys():
            genes_in_network.append(key)

        well_id_in_network_list = []
        for gene in genes_in_network:
            try:
                # gene is prefixed by mtu:, clone_ORF_lookup excludes this prefix
                index = self.clone_ORF_lookup.iloc[:, 1].tolist().index(gene[4:])
                well_id_in_network_list.append(self.clone_ORF_lookup.iloc[index, 0])
            except ValueError:
                pass

        return well_id_in_network_list

    def generate_reaction_drug_effect_dict(self, drug_combo):
        drug_table = self.drug_data.drug_table
        gene_drug_effect_dict = {}
        # generate a dictionary that maps the drug effect to
        # to the well id to gene name (in tb)
        for well_id in self.well_id_in_network_list:
            reverse_index = self.clone_ORF_lookup.iloc[:, 0].tolist().index(well_id)
            gene = self.clone_ORF_lookup.iloc[reverse_index, 1]

            index_of_well_in_drug_table = drug_table.iloc[:, 0].tolist().index(well_id)
            drug_effects = drug_table.iloc[index_of_well_in_drug_table, drug_combo]
            
            summed_drug_effect = sum(drug_effects)
            gene_drug_effect_dict[gene] = summed_drug_effect

        self.gene_drug_effect_dict = gene_drug_effect_dict

        # this section will use the gene to reaction dict
        # and the drug effet dictionary to create a new dictionary
        # which contains reactions as keys and values of effect
        reaction_drug_effect_dict = {}
        for key in self.gene_drug_effect_dict.keys():
            drug_effect = self.gene_drug_effect_dict[key]
            new_key = "mtu:" + key
            reactions = self.genes_to_reactions_dict[new_key]
            for reaction in reactions:
                reaction_drug_effect_dict[reaction] = drug_effect
        
        self.reaction_drug_effect_dict = reaction_drug_effect_dict

        return self.reaction_drug_effect_dict


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
    net = KEGG_Network(pathway_obj)

    lookup = read_cloneID_to_orf_table(path_cloneID_ORF)

    d_u = Determining_Uberedges(drug_data, lookup, net)

    # this section will find the drug effect for each reaction in the list
    drug_combo = [10, 14, 17, 24]

    reac_drug_effect = d_u.generate_reaction_drug_effect_dict(drug_combo)
    
    print(d_u.gene_drug_effect_dict)
    print('###########################################')
    # it may be easier to get stay at the gene level
    # (for integration into other classes)
    print(d_u.reaction_drug_effect_dict)
