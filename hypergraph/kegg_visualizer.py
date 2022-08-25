import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import networkx as nx
from fileIO import read_KGML
from kegg_network import KEGG_Network
import helper_functions
from networkx.drawing.nx_agraph import graphviz_layout


class KEGG_Visualizer:

    """Summary
    
    Attributes
    ----------
    network : TYPE
        Description
    reaction_color_map : TYPE
        Description
    reaction_df : TYPE
        Description
    reaction_graph : TYPE
        Description
    relation_color_map : list
        Description
    relation_graph : TYPE
        Description
    """
    
    def __init__(self, network):
        """Summary
        
        Parameters
        ----------
        network : TYPE
            Description
        """
        self.network = network

        self.relation_df, self.relation_entry_to_name_dict = self.create_relation_df()
        self.relation_name_to_entry_dict = helper_functions.dictionary_reverser(self.relation_entry_to_name_dict)
        
        reac_df, reac_id_dict, sub_id_dict, prod_id_dict = self.create_reaction_df()
        self.reaction_df = reac_df
        self.reaction_id_to_name_dict = reac_id_dict
        self.substrate_id_to_name_dict = sub_id_dict
        self.product_id_to_name_dict = prod_id_dict

        self.relation_graph = self.generate_relation_graph(self.relation_df)
        self.reaction_graph = self.generate_reaction_graph()

        # color  maps will only be set if drug data is added
        self.relation_color_map = None
        self.reaction_color_map = None

    def create_relation_df(self):
        """Summary
        
        Returns
        -------
        TYPE
            Description
        """
        entry_to_name_dict = {}
        entry1_id_list = []
        entry2_id_list = []
        rel_type_list = []
        rel_subtype_list = []

        for relation in self.network.relation_list:
            entry1 = relation.entry1
            entry_to_name_dict[entry1.id] = entry1.name
            entry2 = relation.entry2
            entry_to_name_dict[entry2.id] = entry2.name
            entry1_id_list.append(entry1.id)
            entry2_id_list.append(entry2.id)

            rel_type_list.append(relation.type)
            rel_subtype_list.append(relation.subtypes)

        rel_network_dictionary = {'entry1': entry1_id_list, 'entry2': entry2_id_list, 'rel_type': rel_type_list, 'rel_subtype': rel_subtype_list}
        relation_df = pd.DataFrame(rel_network_dictionary, columns=rel_network_dictionary.keys())

        return relation_df, entry_to_name_dict

    def generate_relation_graph(self, relation_df):
        """Summary
        
        Parameters
        ----------
        relation_df : TYPE
            Description
        
        Returns
        -------
        TYPE
            Description
        """
        G = nx.DiGraph()
        for i in range(0, relation_df.shape[0]):
            edge = relation_df.iloc[i, 0:2].values.tolist()
            label = relation_df.iloc[i, 3]
            G.add_edge(*edge, label=label)

        return G

    def create_reaction_df(self):
        reaction_id_to_name_dict = {}
        substrate_id_to_name_dict = {}
        product_id_to_name_dict = {}
        substrate_list = []
        product_list = []
        reaction_type_list = []  # will be reversible or irreversible
        reaction_id_list = []
        # reaction_subtype_list = []

        for reaction in self.network.reaction_list:
            substrates_entry = reaction.substrates
            products_entry = reaction.products

            temp_subs_list = []
            temp_prods_list = []

            for substrate in substrates_entry:
                temp_subs_list.append(substrate.id)
                substrate_id_to_name_dict[substrate.id] = substrate.name
            for product in products_entry:
                temp_prods_list.append(product.id)
                product_id_to_name_dict[product.id] = product.name

            substrate_list.append(temp_subs_list)
            product_list.append(temp_prods_list)
            reaction_type_list.append(reaction.type)
            reaction_id_to_name_dict[reaction.id] = reaction.name
            reaction_id_list.append(reaction.id)

        react_network_dictionary = {
            'reaction_id': reaction_id_list, 'substrate': substrate_list, 
            'product': product_list, 'reaction_type': reaction_type_list}
        
        reaction_df = pd.DataFrame(
            react_network_dictionary, columns=react_network_dictionary.keys())

        return reaction_df, reaction_id_to_name_dict, substrate_id_to_name_dict, product_id_to_name_dict

    def generate_reaction_graph(self):
        return None

    def plot_graph(self, graph_type):
        """Summary
        
        Parameters
        ----------
        graph_type : TYPE
            Description
        """
        if graph_type == 'relation':
            self.plot_relation_graph()
        elif graph_type == 'reaction':
            self.plot_reaction_graph()
        else:
            print('unrecognized graph type, use relation or reaction')

    def plot_relation_graph(self):
        """Summary
        """

        # setting up the figure (to enable adding a title)
        plt.figure(figsize=(10, 5))
        ax = plt.gca()
        ax.set_title(self.network.name)

        # section for edge weights and colors
        # edge_attributes = nx.get_edge_attributes(self.relation_graph, 'label')  # dictionary
        # colors = []
        # weights = []
        # for edge in self.relation_graph.edges:

        #     # simply increase edge weight:
        #     weights.append(2)

        #     edge_attributes_list = edge_attributes[edge]
        #     attributes = []
        #     # edge attributes are lists of tuples, we unpack them here
        #     for entry in edge_attributes_list:
        #         for thing in entry:
        #             attributes.append(thing)
        #     attributes_str = ' '.join(attributes)
        #     if 'activation' in attributes_str:
        #         # green
        #         colors.append('#a0e2a8')
        #     # elif 'expression' in attributes_str:
        #     #     colors.append('green')
        #     elif 'inhibition' in attributes_str:
        #         # red
        #         colors.append('#d66969')
        #     else:
        #         # balck/grey
        #         colors.append('#383b38')

        # section for node position and size
        my_pos = graphviz_layout(self.relation_graph, prog='neato')
        # my_pos = nx.spring_layout(self.relation_graph, k=0.1, iterations=20, seed=100)
        my_node_size = 350

        label_dict = {}
        for node in self.relation_graph.nodes:
            # for clarity, we only label with the first entry!!!!!!!!!
            # TODO: change this functionality
            
            current_label = self.relation_entry_to_name_dict[node].split(" ")[0]
            if current_label[0:4] == 'hsa:':
                kegg_entrez_id = current_label[4:]
                symbol = helper_functions.entrez_id_to_symbol(self.entrez_ids_df, self.symbol_names_df, int(kegg_entrez_id))
                if symbol == 'None':
                    label_dict[node] = current_label
                else:
                    label_dict[node] = symbol
            else:
                label_dict[node] = current_label

            # label_dict[node] = self.relation_entry_to_name_dict[node].split(" ")[0]

        if self.relation_color_map is None:
            nx.draw(
                self.relation_graph, node_size=my_node_size, pos=my_pos,
                labels=label_dict, with_labels=True, ax=ax)
        else:
            nx.draw(
                self.relation_graph, node_color=self.relation_color_map,
                node_size=my_node_size, pos=my_pos,
                labels=label_dict, with_labels=True, ax=ax)

        # nx.draw_networkx_edge_labels(self.relation_graph, my_pos, font_size=7, edge_labels=nx.get_edge_attributes(self.relation_graph, 'label'), clip_on=False, alpha=0.5)
        # print(nx.get_edge_attributes(self.relation_graph, 'label'))

        plt.show()

    def plot_reaction_graph(self):
        return None

    def add_drug_data(self, drug_data, graph_type):
        """Summary
        
        Parameters
        ----------
        drug_data : dataframe
            should be a single time point of DE drug data
        graph_type : string
            argument should be 'relation' or 'reaction', which graph we want
            to add drug data for
        """
        if graph_type == 'relation':
            annotated_graph = self.annotate_relation_graph(drug_data)
        elif graph_type == 'reaction':
            annotate_graph = self.annotate_reaction_graph(drug_data)
        else:
            print("graph type not recognized, specify relation or reaction")

    def annotate_relation_graph(self, drug_data):
        """Summary
        
        Parameters
        ----------
        drug_data : TYPE
            Description
        
        Returns
        -------
        TYPE
            Description
        """

        # determine which genes are up and which genes are down
        indices_of_network_genes_in_drug_data = helper_functions.extract_indices_from_pathway(drug_data, self.network)

        up_reg_df = drug_data[drug_data.iloc[:, 2] >= 0.1]
        down_reg_df = drug_data[drug_data.iloc[:, 2] <= -0.1]
        # nan_reg_df = drug_data[drug_data.iloc[:, 1].isna()]
        
        regulation_dict = {}
        for entrez_id in up_reg_df.iloc[:, 0]:
            regulation_dict['hsa:' + str(entrez_id)] = 'up'
        for entrez_id in down_reg_df.iloc[:, 0]:
            regulation_dict['hsa:' + str(entrez_id)] = 'down'

        self.relation_color_map = []
        for node in self.relation_graph.nodes:
            node_name = self.relation_entry_to_name_dict[node]
            node_regulation = regulation_dict.get(node_name.split(' ')[0])

            # colors are in hex
            if node_regulation == 'up':
                self.relation_color_map.append('#a0e2a8')
            elif node_regulation == 'down':
                self.relation_color_map.append('#d66969')
            else:
                self.relation_color_map.append('#b5b5b5')

        return None

    def annotate_reaction_graph(self, drug_data):
        return None


if __name__ == '__main__':
    absolutePath = os.path.abspath(__file__)
    fileDirectory = os.path.dirname(absolutePath)
    parentDirectory = os.path.dirname(fileDirectory)
    path_drug_data = os.path.join(
        parentDirectory, "input_files/Multidrug_6hr_Responses_trimmed.xlsx"
    )

    path_KEGG = os.path.join(parentDirectory, "input_files/KEGG_data/")

    file_path = path_KEGG + "mtu00010.xml"

    pathway_obj = read_KGML(file_path)
    network = KEGG_Network(pathway_obj)

    net_vis = KEGG_Visualizer(network)
    net_vis.plot_relation_graph()

    print(len(network.relation_list))
    print(len(network.reaction_list))