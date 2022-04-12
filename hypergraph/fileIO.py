import os
import sys
import pandas as pd


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
    return pd.ExcelFile(path)


def combineSheets(xls):
    """Reads each sheet from given xlsx file and combines them by the first column

    Reads the XLSX file `xls` and creates a list of pandas DataFrames, one for
    each sheet in `xls`. Then concatenates each sheet to the first using the first
    sheets first column.

    Parameters
    ----------
    xls : pd.ExcelFile
        A Pandas ExcelFile that stores the sheets we want to combine

    Returns
    -------
    combinedSheet : pd.DataFrame
        pandas dataframe that contains the information from every sheet in `xls`
    """
    sheets = pd.read_excel(xls,sheet_name=None)
    combinedSheet = pd.concat(sheets, axis=1)
    print(combinedSheet)
    return combinedSheet

    

if __name__ == "__main__":
    absolutePath = os.path.abspath(__file__)
    fileDirectory = os.path.dirname(absolutePath)
    parentDirectory = os.path.dirname(fileDirectory)
    path = os.path.join(
        parentDirectory, "input_files/Multidrug_6hr_Responses_trimmed.xlsx"
    )

    print(readXLSX(path))
    print(path)
    combineSheets(readXLSX(path))
