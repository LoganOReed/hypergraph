import pandas as pd
import os
from fileIO import readXLSX
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
        self.significance_filters = self.create_significance_filters(self.drug_table)
        ##############################

    def create_drug_table(self):
        # this solution is following a stack overflow post
        # https://stackoverflow.com/questions/53935848/how-to-merge-all-data-frames-in-a-dictionary-in-python
        # 'on' determines which col we use a the key for the merge,
        # 'how' is what type of merge, we choose outer
        my_reduce = functools.partial(pd.merge, on='ID_REF', how='outer')
        collated_table = functools.reduce(my_reduce, self.raw_dataframe_dict.values())
        return collated_table 

    def get_drug_list(self):
        # for now, we will treat drugs at different doses as different treatments
        # we just iterate over the keys in the dict and append them to a list
        drug_list = []
        for key in self.raw_dataframe_dict.keys():
            drug_list.append(key)
        return drug_list

    def get_num_replicates(self):
        # this function will take in the dictionary of sheets and
        # return a list that holds the number of replicates for each drug
        num_replicates = []
        for sheet in self.raw_dataframe_dict.values():
            # number of replicates is just the number of columns in each sheet
            # we subtract 1 because the first column is id_ref (ie: gene names)
            num_replicates.append(len(sheet.columns)-1)
        return num_replicates

    def create_average_drug_table(self):
        # create a new dataframe that we will populate with average values
        # we start by setting the first col to the ID_REF of the existing table
        average_drug_table = pd.DataFrame([self.drug_table.ID_REF]).transpose()
        
        for i in range(0, self.num_drugs):
            # build the new column header
            cur_drug = self.drug_list[i]
            cur_col_header = cur_drug + '_avg'

            # get the dataframe for the current drug
            cur_drug_table = self.raw_dataframe_dict[cur_drug]

            # average the values in the columns and update the average_drug_table
            average_drug_table[cur_col_header] = cur_drug_table.iloc[:, 1:].mean(axis=1)

        return average_drug_table

    def create_significance_filters(self, drug_dataframe):
        # drug_dataframe will be a dataframe that has the first col
        # as ID_REF and the rest of the columns will hold gene expression
        # levels

        return 0


if __name__ == '__main__':
    absolutePath = os.path.abspath(__file__)
    fileDirectory = os.path.dirname(absolutePath)
    parentDirectory = os.path.dirname(fileDirectory)
    path = os.path.join(
        parentDirectory, "input_files/Multidrug_6hr_Responses_trimmed.xlsx")

    excel = readXLSX(path)

    drug_data = Drug_Data(excel)

    # print(drug_data.sheet_names)
    # print(drug_data.drug_table)

    # print(type(drug_data.average_drug_table))
    # print(drug_data.average_drug_table.columns)
    # print(drug_data.average_drug_table.shape)

    average_table = drug_data.create_average_drug_table()
    # number of cols and sheet names differ by one, makes sense
    # as the first col is the clone_id
    print(average_table)
    print(len(drug_data.sheet_names))
