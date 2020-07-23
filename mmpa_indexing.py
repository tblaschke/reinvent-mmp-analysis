import pandas as pd
import numpy as np
from rdkit import Chem 
from rdkit.Chem import AllChem
from multiprocessing import Pool, cpu_count
import os
import platform
from glob import glob

OVERWRITE_FILES = False

def applyParallel(df, func):
    df_split = np.array_split(df, cpu_count())
    pool = Pool(cpu_count())
    data = pd.concat(pool.map(func, df_split))
    pool.close()
    pool.join()
    return data

def normalize_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        return Chem.MolToSmiles(mol, isomericSmiles=False)
    else:
        return np.NaN

def normalize_smiles_series(smiles_series):
    return smiles_series.map(normalize_smiles)


rnd_seed = 1234

zcat = "gzcat" if platform.system() == "Darwin" else "zcat"

with open("fragmented_smiles/jobs", 'w') as fd:
    for folder in glob("fragmented_smiles/*_*"):
        foldername = folder.split("/")[-1]
        if os.path.exists(f'MMP/{foldername}/MMP.csv.gz'):
            if not OVERWRITE_FILES:
                print(f"Skipping {foldername} as it seems to already be processed")
                continue          
        os.makedirs(f"MMP/{foldername}", exist_ok=True)
        target = foldername.split("_")[0] 
        files = sorted(glob(f"fragmented_smiles/{foldername}/*.csv.gz")) + sorted(glob(f"fragmented_smiles/{target}/*.csv.gz"))
        files = " ".join(files)
        fd.write(f"{zcat} " + files + " | python mmpa/indexing.py -s -r 0.3333 | gzip > MMP/{}/MMP.csv.gz\n".format(foldername))