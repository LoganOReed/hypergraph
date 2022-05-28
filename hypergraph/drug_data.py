import pandas as pd
import os
from fileIO import readXLSX
import numpy as np
import functools


class Drug_Data:
    def __init__(self, excel_df_dict):
        self.raw_dataframe_dict = excel_df_dict
        self.sheet_names = self.raw_dataframe_dict.keys()

        self.drug_table = self.create_drug_table()
        self.drug_list = self.get_drug_list()
        self.num_drugs = len(self.drug_list)
        self.num_replicates = self.get_num_replicates()
        self.average_drug_table = self.create_average_drug_table()

        # TODO: ############################
        default_sig_level = 1
        self.significance_filters = self.create_significance_filters(
            self.drug_table, default_sig_level
        )
        ##############################

    def create_drug_table(self):
        """This function creates a pd dataframe that holds drug-treated
        expression levels. Data is pulled from self.raw_dataframe_dict
        which is a dictionary where keys are sheet names and values are
        pd dataframes of those sheets

        Returns
        -------
        pd dataframe
            function merges all of the dataframes in self.raw_dataframe_dict
            according to the first column - clone_id
            Format of the final dataframe is: first col is gene clone_ids
            rows correspond to specific gene, cols correspond to drug
            entries are relative gene expression levels (to control)
        """
        # this solution is following a stack overflow post
        # https://stackoverflow.com/questions/53935848/how-to-merge-all-data-frames-in-a-dictionary-in-python
        # 'on' determines which col we use a the key for the merge,
        # 'how' is what type of merge, we choose outer
        my_reduce = functools.partial(pd.merge, on="ID_REF", how="outer")
        collated_table = functools.reduce(my_reduce, self.raw_dataframe_dict.values())
        return collated_table

    def get_drug_list(self):
        """This function simply gets the sheet names from the dictionary

        Returns
        -------
        list
            the return list contains all of the keys from the raw_dataframe_dict
            (the keys are the names of the individual sheets) Each sheet contains
            gene expression levels for treatment with a single drug
        """
        # for now, we will treat drugs at different doses as different treatments
        # we just iterate over the keys in the dict and append them to a list
        drug_list = []
        for key in self.raw_dataframe_dict.keys():
            drug_list.append(key)
        return drug_list

    def get_num_replicates(self):
        """this function takes in a dictionary of dataframes and returns
        the number of columns for each dataframe in the dictionary
        Columns correspond to treatment with the drug (sheet name)
        some experiments (treatments) contain multiple replicates

        Returns
        -------
        list
            list containing the number of replicates for each sheet
        """
        # this function will take in the dictionary of sheets and
        # return a list that holds the number of replicates for each drug
        num_replicates = []
        for sheet in self.raw_dataframe_dict.values():
            # number of replicates is just the number of columns in each sheet
            # we subtract 1 because the first column is id_ref (ie: gene names)
            num_replicates.append(len(sheet.columns) - 1)
        return num_replicates

    def create_average_drug_table(self):
        """the function gives a dataframe containing average drug effect
        we are averaging over the replicates for each drug

        Returns
        -------
        dataframe
            dataframe where rows correspond to gene clone_id and cols
            correspond to drugs. Entries are average (over replicates)
            effect of drug on specific gene
        """
        # create a new dataframe that we will populate with average values
        # we start by setting the first col to the ID_REF of the existing table
        average_drug_table = pd.DataFrame([self.drug_table.ID_REF]).transpose()

        for i in range(0, self.num_drugs):
            # build the new column header
            cur_drug = self.drug_list[i]
            cur_col_header = cur_drug + "_avg"

            # get the dataframe for the current drug
            cur_drug_table = self.raw_dataframe_dict[cur_drug]

            # average the values in the columns and update the average_drug_table
            average_drug_table[cur_col_header] = cur_drug_table.iloc[:, 1:].mean(axis=1)

        return average_drug_table

    def create_significance_filters(self, drug_dataframe, significance_cutoff):
        """this function creates a dataframe where rows correspond to gene id
        and cols correponds to treatment, entries are true/false depending on
        whether gene expression level in drug_dataframe is greater than
        significance_cutoff (in abs values)

        Parameters
        ----------
        drug_dataframe : pd dataframe
            dataframe where first col is clone_id, cols correspond to drugs
            rows correspond to genes, entries are relative gene expression
            levels
        significance_cutoff : double
            value for which we say a drug has a significant effect on a gene

        Returns
        -------
        dataframe
            values of dataframe are True/False - depending on if values are
            above/below the significance level
            first column in clone_ids
        """
        # drug_dataframe will be a dataframe that has the first col
        # as ID_REF and the rest of the columns will hold gene expression
        # levels

        drug_values = drug_dataframe.iloc[:, 1:]  # get everything except the clone id
        clone_ids = drug_dataframe.iloc[:, 0]

        sig_drug_values = np.abs(drug_values) >= significance_cutoff

        # re append the values to the clone ids
        sig_drug_dataframe = pd.concat([clone_ids, sig_drug_values], axis=1)
        return sig_drug_dataframe


if __name__ == "__main__":
    absolutePath = os.path.abspath(__file__)
    fileDirectory = os.path.dirname(absolutePath)
    parentDirectory = os.path.dirname(fileDirectory)
    path = os.path.join(
        parentDirectory, "input_files/Multidrug_6hr_Responses_trimmed.xlsx"
    )

    excel = readXLSX(path)

    drug_data = Drug_Data(excel)

    # print(drug_data.sheet_names)
    # print(drug_data.drug_table)

    # print(type(drug_data.average_drug_table))
    # print(drug_data.average_drug_table.columns)
    # print(drug_data.average_drug_table.shape)

    # print(drug_data.average_drug_table)
    # number of cols and sheet names differ by one, makes sense
    # as the first col is the clone_id
    print(drug_data.significance_filters)
