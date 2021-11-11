import pickle

class PPINetwork:
    def __init__(self, ppi_path, drug_path):
        self.ppi = pickle.load(ppi_path)
        self.drugs = pickle.load(drug_path)

    

    