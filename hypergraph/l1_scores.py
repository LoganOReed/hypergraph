import pandas as pd
import numpy as np
from fileIO import readXLSX, read_KGML, read_cloneID_to_orf_table
from drug_data import Drug_Data
from kegg_network import KEGG_Network
import os


class L1_scores:
    def __init__(self, drug_data):
        self.drug_data = drug_data
        self.drug_table = self.drug_data.drug_table
        self.significance_filters = self.drug_data.significance_filters
        self.drug_list = self.drug_data.drug_list

    def compute_L1(self, drug_combination, network, clone_ORF_lookup):
        """Computes the L1 score for a single drug combination and network
        
        Parameters
        ----------
        drug_combination : list of col indices
            indices corresp to drugs from the drug_table
        network : kegg_network object
            Description
        
        Returns
        -------
        L1 score
            Description
        """

        # convert kegg_network genes to well_ids
        well_id_list = []
        for gene in network.gene_list:
            # gene_lower = gene.lower()
            # table_lower = [x.lower() for x in clone_ORF_lookup.iloc[:, 1]]
            try:
                well_id_index = clone_ORF_lookup.iloc[:, 1].tolist().index(gene[4:])
                well_id_list.append(clone_ORF_lookup.iloc[well_id_index, 0])
            except ValueError:
                # ie: gene isn't in the list, we pass
                # print(gene)
                pass

        # print(well_id_list)

        # apply filter to drug values to ignore insignificant values
        drug_values = self.drug_table.iloc[:, 1:]
        clone_ids = self.drug_table.iloc[:, 0]
        sig_drug_mask = self.significance_filters.iloc[:, 1:]
        sig_drug_values = drug_values[sig_drug_mask]
        sig_drug_table = pd.concat([clone_ids, sig_drug_values], axis=1)

        # find which genes occur the specified pathway
        genes_in_network_index_list = []
        for well_id in well_id_list:
            try:
                gene_in_network_idex = clone_ids.tolist().index(well_id)
                genes_in_network_index_list.append(gene_in_network_idex)
            except ValueError:
                pass
        
        # print(genes_in_network_index_list)

        # compute the sum for the drug combo
        # start with non-abs-val sum over rows (accounts for overlap)
        # finish with column sums in abs val
        sig_drug_combo_values_in_network = sig_drug_values.iloc[genes_in_network_index_list, drug_combination]
        row_sum = sig_drug_combo_values_in_network.sum(axis=1)
        # print(row_sum)
        L1_score = row_sum.abs().sum()

        return L1_score


if __name__ == '__main__':
    absolutePath = os.path.abspath(__file__)
    fileDirectory = os.path.dirname(absolutePath)
    parentDirectory = os.path.dirname(fileDirectory)
    path_drug_data = os.path.join(
        parentDirectory, "input_files/Multidrug_6hr_Responses.xlsx"
    )

    path_KEGG = os.path.join(parentDirectory, "input_files/KEGG_data/")

    path_network = path_KEGG + 'mtu01200.xml'

    path_cloneID_ORF = os.path.join(parentDirectory, 'input_files/clone_to_orf.csv')

    excel = readXLSX(path_drug_data)
    drug_data = Drug_Data(excel)

    pathway_obj = read_KGML(path_network)
    network = KEGG_Network(pathway_obj)

    clone_ORF_lookup = read_cloneID_to_orf_table(path_cloneID_ORF)

    print(clone_ORF_lookup)

    # print(network.gene_list)

    # testing implementation
    drug_combo = [10, 23, 24, 25]
    test = L1_scores(drug_data)
    score = test.compute_L1(drug_combo, network, clone_ORF_lookup)
    print(score)
