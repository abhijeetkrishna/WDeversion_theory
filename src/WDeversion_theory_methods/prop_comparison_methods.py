import os
import urllib.request
import pandas as pd

def load_raw_data(force_download=False):

    # download raw data in folder raw_data if not already present
    webpages = ["https://git.mpi-cbg.de/fuhrmann/wdeversion_paper/-/blob/main/CurvedTM/CurvedTM_notebooks/CellShapeAnalysis/DFallDiscsIncreaselimitcounts.pkl",
                "https://git.mpi-cbg.de/fuhrmann/wdeversion_paper/-/blob/main/CurvedTM/CurvedTM_notebooks/CellShapeAnalysis/DFallDiscslimitcounts.pkl"
                ]

    for webpage in webpages:
        filename = webpage.split("/")[-1]
        if not os.path.isfile("raw_data/"+filename) or force_download:
            os.makedirs("raw_data", exist_ok=True)
            print("Downloading file: {}".format(filename))
            urllib.request.urlretrieve(webpage, "raw_data/"+filename)
        else:
            print("File {} already present".format(filename))  

    # load data
    DFallDiscsIncreaselimitcounts = pd.read_pickle("raw_data/DFallDiscsIncreaselimitcounts.pkl")
    DFallDiscslimitcounts = pd.read_pickle("raw_data/DFallDiscslimitcounts.pkl")

    return(DFallDiscsIncreaselimitcounts, DFallDiscslimitcounts)