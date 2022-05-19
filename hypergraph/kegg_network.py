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
        self.adjacency_matrix = self.create_simple_S_matrix()

    def create_simple_S_matrix(self):
        """creates the S matrix based on the reactions
        from the KGML pathway. Utilizes self.reactions and self.compound_list
        
        We will construct a list of substrates and products that are in
        parallel - where the same index will denote a reaction/edge

        We will then iterate through those lists, constructing the simple edge
        S matrix. Substrates will get a value of -1, while products get a value
        of +1.

        The S matrix will be built one column at a time (cols in S matrix
        correspond to reactions/edges)

        Returns
        -------
        2d matrix
            
        """

        # we iterate throgh the reactions, building product list
        # and substrate list, note, if reaction is reversible then we add
        # substrates to the product list and products to the substrates list
        # as well (ie: the opposite direction)
        substrate_list = []
        product_list = []
        for reaction in self.reactions:
            substrates = reaction.substrates
            products = reaction.products
            
            substrate_list.append(substrates)
            product_list.append(products)

            if reaction.type == 'reversible':
                product_list.append(substrates)
                substrate_list.append(products)    

        # pre-processing on the substrates and products
        # note: these lists are after duplication for reversible reactions
        # sub_name_expanded and prod_name_expanded are in parallel - corresp to
        # reactions
        sub_name_expanded_list = []
        for substrate in substrate_list:
            sub_name_list = [sub._names for sub in substrate]  # extract the names
            temp_list = []
            for element in sub_name_list:
                temp_list += element
            sub_name_expanded_list.append(temp_list)

        prod_name_expanded_list = []
        for product in product_list:
            prod_name_list = [prod._names for prod in product]
            temp_list = []
            for element in prod_name_list:
                temp_list += element
            prod_name_expanded_list.append(temp_list)

        # create the simple edge S matrix:
        num_edges = len(sub_name_expanded_list)
        num_compounds = len(self.compound_list)
        simple_S_matrix = np.zeros((num_compounds, 1))  # first col of zeros for shape
        for i in range(0, num_edges):
            # get the subs/prods for the reaction
            substrates = sub_name_expanded_list[i]
            products = prod_name_expanded_list[i]

            # initialize empty col array
            col = np.zeros((num_compounds, 1))

            # find the index of the substrates in the compound list
            indices_of_subs = np.where(np.isin(self.compound_list, substrates))[0]
            indices_of_prods = np.where(np.isin(self.compound_list, products))[0]

            col[indices_of_subs] = -1
            col[indices_of_prods] = 1

            simple_S_matrix = np.hstack((simple_S_matrix, col))
        simple_S_matrix = simple_S_matrix[:, 1:]  # get rid of col of zeros

        print(simple_S_matrix.shape)

        # sample code for demonstrating odd behavior
        for reaction in self.reactions:
            test = [entry._names for entry in reaction.substrates]
            print('##########################################')
            print(reaction)
            print('reaction.substrates:', reaction.substrates)
            print('substrates._names: ')
            for thing in reaction.substrates:
                print(thing._names)

        return simple_S_matrix

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

    # print('genes:')
    # print(network.gene_list)
    # print('#######################')
    # print('compounds:')
    # print(network.compound_list)
