import pandas as pd
import numpy as np
import os
from kegg_network import KEGG_Network
from kegg_visualizer import KEGG_Visualizer
from fileIO import readXLSX, read_KGML, read_cloneID_to_orf_table
import helper_functions
import networkx as nx
import matplotlib.pyplot as plt


class Network_Completion:

    """Summary
    
    Attributes
    ----------
    edge_df : pandas DataFrame
        dataframe that contains reactions/edges
        includes: reac_id, subs, prods, type
    edge_df_columns : pandas series
        column names for the dataframe
    network : KEGG_Network
        network that we wish to complete
    partial_complete_network : pandas DataFrame
        edge dataframe with extra edges
    """
    
    def __init__(self, network, edge_df):
        """Initializer for the Network_Completion class
        
        Parameters
        ----------
        network : KEGG_Network
            network to be completed
        edge_df : pandas DataFrame
            dataframe that contains information about all edges
            within the network
        """
        self.network = network
        self.edge_df = edge_df
        self.edge_df_columns = self.edge_df.columns
        self.partial_complete_network = self.partial_complete_network()

    def partial_complete_network(self):
        """Partially completes the network
        
        Returns
        -------
        pandas DataFrame
            dataframe that contains extra required edges for the sources
            and sinks into the network
        """
        # get unique substrates
        substrate_set = set()
        for entry in self.edge_df.substrate:
            for sub in entry:
                substrate_set.add(sub)

        # get unique products
        product_set = set()
        for entry in self.edge_df['product']:
            for prod in entry:
                product_set.add(prod)

        # get a list of all compounds in the network
        all_compounds = substrate_set.union(product_set)

        # if a substrate doesn't appear in products
        # then it never has an incoming edge, so it needs a source
        nodes_that_need_source = []
        for sub in substrate_set:
            if sub not in product_set:
                nodes_that_need_source.append([sub])

        # construct a dictionary to append to the edge list
        # for 'reactions' that correspond to sources
        # we asssume column names here
        temp_products = nodes_that_need_source  # these substrates need sources
        temp_substrates = [[-1]] * len(nodes_that_need_source)  # [-1] will correspond to the source
        temp_reaction_id = [-1] * len(nodes_that_need_source)
        temp_reaction_type = ['source'] * len(nodes_that_need_source)
        source_dict = {
            'reaction_id': temp_reaction_id, 'substrate': temp_substrates,
            'product': temp_products, 'reaction_type': temp_reaction_type}

        source_df = pd.DataFrame(source_dict)

        # if a product doesn't appear in the substrate
        # then it never has an outgoing edge, so it needs a sink
        nodes_that_need_sink = []
        for prod in product_set:
            if prod not in substrate_set:
                nodes_that_need_sink.append([prod])

        # construct a dictionary to append to the edge list
        # for 'reactions' that correspond to sinks
        # we asssume column names here
        temp_products = [[-2]] * len(nodes_that_need_sink)  # [-2] will be out sink
        temp_substrates = nodes_that_need_sink  # these nodes need a sink
        temp_reaction_id = [-2] * len(nodes_that_need_sink)
        temp_reaction_type = ['sink'] * len(nodes_that_need_sink)
        sink_dict = {
            'reaction_id': temp_reaction_id, 'substrate': temp_substrates,
            'product': temp_products, 'reaction_type': temp_reaction_type}

        sink_df = pd.DataFrame(sink_dict)

        # construct partially completed network
        partially_completed_network = pd.concat([self.edge_df, source_df, sink_df], axis=0)
        return partially_completed_network


if __name__ == '__main__':
    absolutePath = os.path.abspath(__file__)
    fileDirectory = os.path.dirname(absolutePath)
    parentDirectory = os.path.dirname(fileDirectory)
    path_drug_data = os.path.join(
        parentDirectory, "input_files/Multidrug_6hr_Responses.xlsx"
    )

    path_KEGG = os.path.join(parentDirectory, "input_files/KEGG_data/")
    path_cloneID_ORF = os.path.join(parentDirectory, "input_files/clone_to_orf.csv")
    cloneID_ORF = read_cloneID_to_orf_table(path_cloneID_ORF)

    clone_id_to_gene = {}
    for i in range(0, cloneID_ORF.shape[0]):
        clone = cloneID_ORF.iloc[i, 0]
        gene = cloneID_ORF.iloc[i, 1]
        
        clone_id_to_gene[clone] = gene

    file_path = path_KEGG + "mtu00010.xml"

    # we create a pathway object, then create a network object
    # Finally, make a KEGG_Visualizer object based on the network
    pathway_obj = read_KGML(file_path)
    network = KEGG_Network(pathway_obj)

    genes_to_reactions_dict = helper_functions.genes_to_reactions(network)
    
    # section for creating dictionaries that map
    # kegg id to gene names
    genes_to_node_id = {}
    node_id_to_gene = {}
    for entry in network.reaction_entries:
        genes_to_node_id[entry.name] = entry.id
        for name in entry.name.split(' '):
            node_id_to_gene[entry.id] = name[4:]

    net_vis = KEGG_Visualizer(network, clone_id_to_gene)
    edge_obj = net_vis.reaction_df

    net_comp = Network_Completion(network, edge_obj)
    # expand the hyper edges wihtin the network
    expanded_df = net_vis.expand_reaction_df(net_comp.partial_complete_network)
    # create a networkx graph based on expanded network
    G = net_vis.generate_reaction_graph(expanded_df)
