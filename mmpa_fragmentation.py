import pandas as pd
import numpy as np
from rdkit import Chem 
from rdkit.Chem import AllChem
from multiprocessing import Pool, cpu_count
import os
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

outfiles = []
for file in glob("to_fragment/*/*.smi"):
    
    foldername = file.split("/")[1]
    filename = file.split("/")[-1]
    filename = filename.split(".")[0]
    
    
    splitLen = 1000 
    outputfolder = f"to_fragment/splits/{foldername}"
    os.makedirs(outputfolder, exist_ok=True)
    os.makedirs(f"fragmented_smiles/{foldername}/", exist_ok=True)
    outputBase = filename 

    if os.path.exists(f'fragmented_smiles/{foldername}/{outputBase.replace("_to_fragment","_fragmented")}.part.1.csv.gz'):
        if not OVERWRITE_FILES:
            print(f"Skipping {foldername} as it seems to already be processed")
            continue
        
    # This is shorthand and not friendly with memory, but it works.
    input = open(file, 'r').read().split('\n')

    at = 1
    for lines in range(0, len(input), splitLen):
        outputData = input[lines:lines+splitLen]        
        outfiles.append(foldername + "/" + outputBase + ".part." + str(at) + '.smi')
        output = open(outputfolder + "/" + outputBase + ".part." + str(at) + '.smi', 'w')
        output.write('\n'.join(outputData))
        output.close()

        # Increment the counter
        at += 1

with open("to_fragment/jobs", "w") as fd: #for parallel to fragment on all cores
    for outfile in outfiles:
        fd.write("python mmpa/rfrag.py < to_fragment/splits/{} | gzip > fragmented_smiles/{}\n".format(outfile, outfile.replace("_to_fragment","_fragmented").replace(".smi",".csv.gz")))