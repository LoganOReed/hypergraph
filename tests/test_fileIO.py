import os
import pandas as pd


def test_filePath():
    """Tests file path"""
    # Replace correctPath with your actual path
    # I don't know if the os changes how this code functions
    correctPath = "/home/occam/Documents/code/hypergraph/input_files/Multidrug_6hr_Responses_trimmed.xlsx"
    absolutePath = os.path.abspath(__file__)
    fileDirectory = os.path.dirname(absolutePath)
    parentDirectory = os.path.dirname(fileDirectory)
    path = os.path.join(
        parentDirectory, "input_files/Multidrug_6hr_Responses_trimmed.xlsx"
    )
    assert path == correctPath


def test_fileIO():
    """Tests fileIO.py being able to open a file"""

    # Sheet names for trimmed 6hr responses
    sheets = [
        "PZA_6hr_0.12mg_mL",
        "PZA_6hr_1.2mg_mL",
        "EMB_6hr_10ug_mL",
        "Rifampicin_6hr_0.2ug_mL",
    ]

    absolutePath = os.path.abspath(__file__)
    fileDirectory = os.path.dirname(absolutePath)
    parentDirectory = os.path.dirname(fileDirectory)
    path = os.path.join(
        parentDirectory, "input_files/Multidrug_6hr_Responses_trimmed.xlsx"
    )
    # check that pandas opened a file
    xls = pd.ExcelFile(path)
    assert isinstance(xls, pd.ExcelFile)

    # check that it can read in sheets
    for s in sheets:
        data = pd.read_excel(xls, sheet_name=s)
        print(data)
        assert isinstance(data, pd.DataFrame)
