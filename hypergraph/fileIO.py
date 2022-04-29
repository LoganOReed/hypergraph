import os
import sys
import pandas as pd
from Bio.KEGG.REST import kegg_info, kegg_list, kegg_get
from Bio.KEGG.KGML import KGML_parser
from Bio.Graphics.KGML_vis import KGMLCanvas


def readXLSX(path):
    """Reads XLSX file into Pandas ExcelFile

    Reads the XLSX file at `path` and returns a pandas ExcelFile

    Parameters
    ----------
    path : os.path
        A path to xlsx file

    Returns
    -------
    out : pd.ExcelFile
        A pandas ExcelFile with the same data as `path`
    """
    return pd.read_excel(path, None)


def combineSheets(xls):
    """Reads each sheet from given xlsx file and combines them by the first column

    Reads the XLSX file `xls` and creates a list of pandas DataFrames, one for
    each sheet in `xls`. Then concatenates each sheet to the first using the first
    sheets first column.

    Parameters
    ----------
    xls : dictionary
        A dictionary where keys are sheetnames and values are pd.Dataframes of sheets

    Returns
    -------
    combinedSheet : pd.DataFrame
        pandas dataframe that contains the information from every sheet in `xls`
    """
    sheets = xls.values()
    combinedSheet = pd.concat(sheets, axis=1)
    print(combinedSheet)
    return combinedSheet


def draw_kegg_map(path, map_id):
    """ Render a local PDF of a KEGG map with the passed map ID
    """
    # Get the background image first
    pathway = KGML_parser.read(kegg_get(map_id, "kgml"))
    canvas = KGMLCanvas(pathway, import_imagemap=True)
    img_filename = "%s.pdf" % map_id
    canvas.draw(path + img_filename)


def fetch_KGML_file(path, map_id):
    """this function uses biopython modules to fetch the kgml file
    associated with map_id and saves that kgml file as an xml file
    in the specified path
    
    Parameters
    ----------
    path : string
        path to the KEGG_data folder
    map_id : string
        KEGG pathway identifier, eg: mtu01200 for central carbon metabolism
    """
    pathway = kegg_get(map_id, 'kgml')
    kgml_text = pathway.read()
    file_name = map_id + ".xml"
    with open(path+file_name, 'w') as file:
        file.write(kgml_text)
    file.close()


def get_possible_pathways(organism):
    """Helper function, if we want to see possible pathways for
    a specific organism, this function will display them along
    with a short description
    
    Parameters
    ----------
    organism : string
        KEGG identifier for an organism. mtu is mycobacterium tuberculosis
        hsa is homo sapiens (human)
    
    Returns
    -------
    string
        returns a string that holds the KEGG identifiers for all pathways
        for the specified organism
    """
    pathways_text = kegg_list('pathway', organism)
    pathways = pathways_text.read()

    return pathways


def read_KGML(path):
    """reads .xml file into a dictionary
    
    Parameters
    ----------
    path : string
        path to the .xml file
    
    Returns
    -------
    dicitonary
        python dictionary where 
    """
    KEGG_info = KGML_parser.read(open(path, 'r'))
    return KEGG_info


if __name__ == "__main__":
    absolutePath = os.path.abspath(__file__)
    fileDirectory = os.path.dirname(absolutePath)
    parentDirectory = os.path.dirname(fileDirectory)
    path_drug_data = os.path.join(
        parentDirectory, "input_files/Multidrug_6hr_Responses_trimmed.xlsx"
    )

    path_KEGG = os.path.join(parentDirectory, "input_files/KEGG_data/")

    # print(readXLSX(path))
    # print(path)
    # combineSheets(readXLSX(path))

    pathway_list = ['mtu01200', 'mtu00010']
    for pathway in pathway_list:
        print('skipping over the loading')
        # fetch_KGML_file(path_KEGG, pathway)
