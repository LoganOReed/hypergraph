from kegg_network import KEGG_Network
from fileIO import read_KGML, read_cloneID_to_orf_table, readXLSX
from kegg_visualizer import KEGG_Visualizer
from drug_data import Drug_Data
import helper_functions
import os


class Network_Integration:

    """Class for combining reaction networks organism networks
    
    Attributes
    ----------
    combined_network : KEGG_Network
        baseline network with genes (according to passed organism network)
    org_network : KEGG_Network
        organism-specific network, has genes
    rn_network : KEGG_Network
        baseline network, no genes
    """
    
    def __init__(self, org_network, rn_network):
        self.org_network = org_network
        self.rn_network = rn_network

        rn_network_with_genes = rn_network
        rn_network_with_genes.genes = org_network.genes
        rn_network_with_genes.gene_list = org_network.generate_gene_list

        self.combined_network = rn_network_with_genes

    def get_combined_network(self):
        """getter for combined network
        
        Returns
        -------
        KEGG_Network
            baseline reaction network with genes
        """
        return self.combined_network


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

    # only grab the first ten rows of the average table
    drug_data = all_drug_data.average_drug_table.iloc[:, [0, 10]]

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

    mtu_net_vis = KEGG_Visualizer(mtu_network, clone_id_to_gene)
    mtu_net_vis.add_drug_data(drug_data, graph_type='reaction')
    mtu_net_vis.plot_graph(graph_type='reaction')

    combined_net_vis = KEGG_Visualizer(combined_network, clone_id_to_gene)
    combined_net_vis.add_drug_data(drug_data, graph_type='reaction')
    combined_net_vis.plot_graph(graph_type='reaction')
