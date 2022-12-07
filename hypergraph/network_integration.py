from kegg_network import KEGG_Network
from fileIO import read_KGML
from kegg_visualizer import KEGG_visualizer
import os

if __name__ == '__main__':
    absolutePath = os.path.abspath(__file__)
    fileDirectory = os.path.dirname(absolutePath)
    parentDirectory = os.path.dirname(fileDirectory)
    path_drug_data = os.path.join(
        parentDirectory, "input_files/Multidrug_6hr_Responses_trimmed.xlsx"
    )

    path_KEGG = os.path.join(parentDirectory, "input_files/KEGG_data/")

    mtu_file_path = path_KEGG + "mtu01200.xml"
    rn_file_path = path_KEGG + "rn01200.xml"

    mtu_pathway_obj = read_KGML(mtu_file_path)
    mtu_network = KEGG_Network(mtu_pathway_obj)

    rn_pathway_obj = read_KGML(rn_file_path)
    rn_network = KEGG_Network(rn_pathway_obj)

    print(mtu_network.reaction_list)
    print(rn_network.reaction_list)

