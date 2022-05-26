import pandas as pd
import numpy as np
from fileIO import readXLSX, read_KGML, read_cloneID_to_orf_table
from drug_data import Drug_Data
from kegg_network import KEGG_Network
import os


def compute_L1(drug_table, network_list):

    return 0


if __name__ == '__main__':
    absolutePath = os.path.abspath(__file__)
    fileDirectory = os.path.dirname(absolutePath)
    parentDirectory = os.path.dirname(fileDirectory)
    path_drug_data = os.path.join(
        parentDirectory, "input_files/Multidrug_6hr_Responses_trimmed.xlsx"
    )

    path_KEGG = os.path.join(parentDirectory, "input_files/KEGG_data/")

    path_network = path_KEGG + 'mtu01200.xml'

    path_cloneID_ORF = os.path.join(parentDirectory, 'input_files/clone_to_orf.csv')

    excel = readXLSX(path_drug_data)
    drug_data = Drug_Data(excel)

    pathway_obj = read_KGML(path_network)
    network = KEGG_Network(pathway_obj)

    clone_ORF_lookup = read_cloneID_to_orf_table(path_cloneID_ORF)
