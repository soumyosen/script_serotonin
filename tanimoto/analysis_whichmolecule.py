# run the script though Schrodinger's wrapper to fix the import part
# use -p to prepare the submit script in a new output_folder
#  $SCHRODINGER/run ./run_zinc_tanimoto.py -p -f output_folder -r /sc/arion/projects/H_filizm02a/davide/tanimoto_zinc/reference.fp

# run the tanimoto calculations in the output_folder
# sh ./jobs.sh

#  analyse the results
#  $SCHRODINGER/run ./run_zinc_tanimoto.py -a -f output_folder -r /sc/arion/projects/H_filizm02a/davide/tanimoto_zinc/reference.fp


import os
import warnings
import sys, getopt
import subprocess
import re
import pandas as pd
import numpy as np

from schrodinger import structure
from schrodinger.structutils import analyze
from schrodinger import adapter

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw


def get_tanimoto_files(outputfolder):
     '''
     get the list of tanimoto files
     '''
     topd, dirs, files = next(os.walk(f'{outputfolder}'))
     tfiles = [f'{filename}' for filename in files if re.match(r"[A-Z][A-Z].tnm", filename)]
     ret = [(re.sub(".tnm","",filename), f'{outputfolder}/{filename}') for filename in sorted(tfiles)]
     print(f'{len(ret)} tanimoto files ')
     return(ret)

def read_tanimotos(filelist):
     alldata = []
     for tr, filen in filelist:
         df = pd.read_csv(filen)
         df['tranche']=tr
         df['index']=np.arange(0,df.shape[0])
         alldata.append(df)
     alld = pd.concat(alldata)
     return alld


if __name__ == "__main__":

    #read tanimoto data
    tfiles = get_tanimoto_files("tc_output")
    dff = read_tanimotos(tfiles)
    print(f'all tanimotos read {dff.shape}')
    print(dff.head())
    list_of_mols = ["ZINC000000001137", "ZINC000031483558", "ZINC000000601293", "ZINC000029489118"]
    for i in list_of_mols:
        dff1 = dff[dff["canvassim"]==i]
        print(dff1)
