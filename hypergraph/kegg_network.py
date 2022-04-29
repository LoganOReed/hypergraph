import numpy as np
from fileIO import read_KGML
import os


class KEGG_Network:

    def __init__(self, pathway_object):
        self.pathway = pathway_object

        # different types of edges
        self.reactions = self.pathway.reactions  # iterable of KGML_pathway.Reaction
        self.reaction_entries = self.pathway.reaction_entries  # seems to relate to images
        self.relations = self.pathway.relations  # seems to be empty?

        # different types of nodes
        self.genes = self.pathway.genes  # unsure of difference to compounds
        self.gene_list = self.generate_gene_list()

        self.compounds = self.pathway.compounds
        self.compound_list = self.generate_compound_list()

        self.maps = self.pathway.maps

        # TODO:
        self.adjacency_matrix = self.create_adjacency_matrix()

    def create_adjacency_matrix(self):
        """creates an adjacency matrix based on the reactions
        from the KGML pathway. Utilizes self.reactions
        
        Returns
        -------
        TYPE
            Description
        """
        return 0

    def generate_gene_list(self):
        """generates a list of genes that appear in the network
        uses self.genes (a Bio.KEGG object) to generate the list
        
        Returns
        -------
        list
            a list containing all of the genes that appear in the network
        """
        gene_list = []

        for entry in self.genes:
            names = entry.name.split(' ')
            for name in names:
                gene_list.append(name)
        return gene_list

    def generate_compound_list(self):
        """generates a list of the compounds that appear int he network
        uses self.compounds (a Bio.KEGG object) to generate the list
        
        Returns
        -------
        list
            a list containing all of the compounds that appear in the network
        """
        compound_list = []

        for entry in self.compounds:
            names = entry.name.split(' ')
            for name in names:
                compound_list.append(name)
        return compound_list


if __name__ == '__main__':

    absolutePath = os.path.abspath(__file__)
    fileDirectory = os.path.dirname(absolutePath)
    parentDirectory = os.path.dirname(fileDirectory)
    path_drug_data = os.path.join(
        parentDirectory, "input_files/Multidrug_6hr_Responses_trimmed.xlsx"
    )

    path_KEGG = os.path.join(parentDirectory, "input_files/KEGG_data/")

    file_path = path_KEGG + 'mtu01200.xml'
    
    pathway_obj = read_KGML(file_path)
    network = KEGG_Network(pathway_obj)

    # information about these objects:
    # https://biopython.org/docs/1.76/api/Bio.KEGG.KGML.KGML_pathway.html
    # print(network.pathway)
    # for entry in network.reactions:
    #     print(type(entry))
    print('genes:')
    print(network.gene_list)
    print('#######################')
    print('compounds:')
    print(network.compound_list)