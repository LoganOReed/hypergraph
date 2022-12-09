from kegg_network import KEGG_Network
from fileIO import read_KGML, read_cloneID_to_orf_table, readXLSX
from kegg_visualizer import KEGG_Visualizer
from network_integration import Network_Integration
from network_completion import Network_Completion
from drug_data import Drug_Data
import helper_functions
import os
from determining_uberedges import Determining_Uberedges
import numpy as np


if __name__ == '__main__':
    absolutePath = os.path.abspath(__file__)
    fileDirectory = os.path.dirname(absolutePath)
    parentDirectory = os.path.dirname(fileDirectory)
    path_drug_data = os.path.join(
        parentDirectory, "input_files/Multidrug_6hr_Responses.xlsx"
    )

    # import the drug data into a Drug_Data object
    excel = readXLSX(path_drug_data)
    all_drug_data = Drug_Data(excel)
    averaged_drug_data = all_drug_data.average_drug_table.iloc[:, [0, 10]]

    # only grab the first ten rows of the average table
    # drug_data = all_drug_data.average_drug_table.iloc[:, [0, 10]]

    path_KEGG = os.path.join(parentDirectory, "input_files/KEGG_data/")
    path_cloneID_ORF = os.path.join(parentDirectory, "input_files/clone_to_orf.csv")
    cloneID_ORF = read_cloneID_to_orf_table(path_cloneID_ORF)

    # map form GEO clone id to the gene names
    clone_id_to_gene = {}
    for i in range(0, cloneID_ORF.shape[0]):
        clone = cloneID_ORF.iloc[i, 0]
        gene = cloneID_ORF.iloc[i, 1]
        
        clone_id_to_gene[clone] = gene

    mtu_file_path = path_KEGG + "mtu01200.xml"
    rn_file_path = path_KEGG + "rn01200.xml"

    # create the network objects
    mtu_pathway_obj = read_KGML(mtu_file_path)
    mtu_network = KEGG_Network(mtu_pathway_obj)

    rn_pathway_obj = read_KGML(rn_file_path)
    rn_network = KEGG_Network(rn_pathway_obj)

    # section for creating dictionaries that map
    # kegg id to gene names
    mtu_genes_to_node_id = {}
    mtu_node_id_to_gene = {}
    for entry in mtu_network.reaction_entries:
        mtu_genes_to_node_id[entry.name] = entry.id
        for name in entry.name.split(' '):
            mtu_node_id_to_gene[entry.id] = name[4:]

    # map from the genes to the reactions
    genes_to_reactions_dict = helper_functions.genes_to_reactions(mtu_network)

    # create the combined network
    net_integration = Network_Integration(mtu_network, rn_network)
    combined_network = net_integration.combined_network

    # create the visualizer for the combined network
    combined_net_vis = KEGG_Visualizer(combined_network, clone_id_to_gene)
    # combined_net_vis.add_drug_data(averaged_drug_data, graph_type='reaction')
    # combined_net_vis.plot_graph(graph_type='reaction')

    # ge the reaction df from the combined network visualizer
    combined_reaction_df = combined_net_vis.reaction_df

    # complete the network based on the combined network
    # we use the reaction_df from the visualizer, should be refactored at some point
    net_comp = Network_Completion(combined_network, combined_reaction_df)

    # get the completed reaction df from the Network_Completion object
    simple_completed_df = net_comp.partial_complete_network
    
    # determine the uberedges in the combined network based on the drug data
    determine_net_uber_edges = Determining_Uberedges(all_drug_data, cloneID_ORF, combined_network)

    # this section will find the drug effect for each reaction in the list
    drug_combo = [10, 14, 17, 24]

    reac_drug_effect = determine_net_uber_edges.generate_reaction_drug_effect_dict(drug_combo)

    # print(reac_drug_effect)
    # print(simple_completed_df)

    # create a mapping between the reaction name and the reaction id
    # there may be an issue here, we are losing in our mapping
    # each reaction can have more than one name, but only 1 id
    reaction_name_to_id = {}
    for reaction in combined_network.reactions:
        names = reaction.name
        reac_id = reaction.id
        for name in names.split(" "):
            reaction_name_to_id[name] = reac_id

    # print(reaction_name_to_id)

    # create a mapping between the reaction id and the drug effect
    # there may be an issue here, we a re losing in our mapping
    reaction_id_to_drug_effect = {}
    for key in reac_drug_effect.keys():
        drug_effect = reac_drug_effect[key]
        # reac_id = reaction_name_to_id[key]
        reac_id = reaction_name_to_id.get(key, 0)

        if reac_id == 0:
            pass
        else:
            reaction_id_to_drug_effect[reac_id] = drug_effect

    # print(reaction_id_to_drug_effect)

    # we loose a couple....
    print(len(reac_drug_effect.keys()))  # 78
    print(len(reaction_id_to_drug_effect.keys()))  # 52

    # add a row to the dataframe to label if a reaction is affected
    # by an uberedge
    is_uber_edge = []
    # uber_edge_vals = []
    for i in range(0, simple_completed_df.shape[0]):
        reaction_id = simple_completed_df.iloc[i, 0]
        if reaction_id in reaction_id_to_drug_effect.keys():
            is_uber_edge.append(1)
            # uber_edge_vals.append(reaction_id_to_drug_effect[reaction_id])
        else:
            is_uber_edge.append(0)
            # uber_edge_vals.append(np.nan)

    complete_with_uber_df = simple_completed_df
    complete_with_uber_df["is_uber"] = is_uber_edge

    print(complete_with_uber_df)
        
