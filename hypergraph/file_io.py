import os
import pandas as pd


def readXLSX(path):
    """Reads XLSX file into Pandas DataFrame

    Reads the XLSX file at `path` and returns a pandas DataFrame

    Parameters
    ----------
    path : os.path
        A path to xlsx file

    Returns
    -------
    out : pd.DataFrame
        A pandas dataframe with the same data as `path`
    """
    # specifying None returns a dictionary where keys are sheet names
    # and values are dataframes of the corresponding sheet
    return pd.read_excel(path, None)


if __name__ == "__main__":
    absolutePath = os.path.abspath(__file__)
    fileDirectory = os.path.dirname(absolutePath)
    parentDirectory = os.path.dirname(fileDirectory)
    path = os.path.join(
        parentDirectory, "input_files/Multidrug_6hr_Responses_trimmed.xlsx"
    )
    print(readXLSX(path))
