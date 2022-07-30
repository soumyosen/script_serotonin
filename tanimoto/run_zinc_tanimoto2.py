
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

def get_zincfp_tranches(folder):
     '''
     get the list of zinc fingerprint files 
     '''
     topd, dirs, files = next(os.walk(f'{folder}'))
     fpfiles = [f'{filename}' for filename in files if re.match(r"[A-Z][A-Z]_out.fp", filename)]
     ret = [(re.sub("_out.fp","",filename), f'{folder}/{filename}') for filename in sorted(fpfiles)]
     print(f'{len(ret)} fp files ')
     return(ret)

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


dir0 = "/sc/arion/projects/H_filizm02a/ZINC20_3d/fingerprints/20210406"
canvascommand = '"$SCHRODINGER/utilities/canvasFPMatrix" -metric tanimoto'
jobparams = '-HOST Chimera-Long-lic03-2020-42:32'
tmpparams = '-TMPDIR /sc/arion/scratch/sens05/.tmp_schrodinger/'

def  write_jobs(fpzinc, outputfolder, referencefile, params=''):
     cmd1 = [f'-ifp {filen} -ifp2 {referencefile} -ocsv {tr}.tnm -JOB tanimoto_{tr}' for tr,filen in fpzinc]
     commands = [canvascommand + ' ' + params + ' ' + cmd + ' ' + jobparams + ' ' + tmpparams for cmd in cmd1]
     return commands

def extract_ligs(ligand_names, filename):
      '''
      extracts the ligands based on name
      '''
      extracted_list = []
      with structure.StructureReader(filename) as reader:
           for i, st in enumerate(reader):
                #print(st.title, st.unique_smiles) 
                if st.title in ligand_names:
                       print(i)
                       smi = analyze.generate_smiles(st)
                       print(st.title, smi)
                       extracted_list.append(st)
      return extracted_list


def extract_ligs_index(ligand_index, filename, properties=None):
      '''
      extracts the ligands based on name
      '''
      extracted_list = []
      with structure.StructureReader(filename) as reader:
           for i, st in enumerate(reader):
                #print(st.title, st.unique_smiles) 
                if i in ligand_index:
                       #print(i, st.title)
                       smi = analyze.generate_smiles(st)
                       print(i, st.title, smi)
                       if properties:
                            k = ligand_index.index(i)
                            for prop in properties:
                                  st.property[prop] = properties[prop][k]
                                  #print(st.property[prop])
                       extracted_list.append(st)
      return extracted_list





if __name__ == "__main__":

   prep_submit = False
   analyse=False
   outputfolder=None
   referencefile=None
   params=''

   try:
        opts, args = getopt.getopt( sys.argv[1:],"apf:r:s:", ["analyse","prepare","folder=","reference=","settings="],)
   except getopt.GetoptError:
        print('prep.py -p -f <foldername> -r <reference>')
   #print(opts)
   for opt, arg in opts:
        if opt in ["-a","--analyse"]:
              analyse = True
        if opt in ["-p","--prepare"]:
              prep_submit= True
        if opt in ["-f","--folder"]:
              outputfolder=arg
        if opt in ["-r","--reference"]:
              referencefile=arg
        if opt in ["-s","--settings"]:
              params=arg

   if prep_submit:

       fpzinc = get_zincfp_tranches(dir0)
       print(f'job will be prepared in folder {outputfolder}')
       print(f'reference fp is {referencefile}')

       if len(fpzinc)>0:
          jobcommands = write_jobs(fpzinc, outputfolder, referencefile, params)
          if not os.path.exists(f'{outputfolder}'):
              os.makedirs(f'{outputfolder}')
              with open(f'{outputfolder}/jobs.sh', 'w') as fileh:
                   fileh.writelines("%s\n" % command for command in jobcommands)
              print(f'command file saved in {outputfolder}/jobs.sh')
          else:
              warnings.warn("output folder exist!")
       else:
            warnings.warn("No files to process", RuntimeWarning)
   

   if analyse:

      lower_cutoff=0.3
      ligpreppdir = '/sc/arion/projects/H_filizm02a/ZINC20_3d//20210406_glide_ready_zinc//'

      #read tanimoto data
      tfiles = get_tanimoto_files(outputfolder)
      dff = read_tanimotos(tfiles)
      print(f'all tanimotos read {dff.shape}')

      #extract the names of the reference ligands
      refs = [col for col in dff.columns if not col in ['canvassim','tranche','index']]
      print(f'reference ligands {refs}')
    
      # sort by the similarity to the first ref ligands
      dff1 =  dff.sort_values(by=refs[0], ascending=False)
      #print(dff1.head())

      #print distribution
      dff1['tcbin'] = pd.cut(dff1[refs[0]], bins = np.arange(0,51)/50)
      with pd.option_context('display.max_rows', None, 'display.max_columns', None):  
           print( dff1['tcbin'].value_counts(sort=False))

      #### modified by Soumyo to fix the removal of few structures of a ligand due to hard TC cutoff 
      dff2 = dff1[dff1[refs[0]] > lower_cutoff]
      dff2_lig_list = dff2["canvassim"].unique().tolist()
      df_extract = dff1[dff1["canvassim"].isin(dff2_lig_list)]
      #df_extract = dff1[dff1[refs[0]] > lower_cutoff]
      #############################################################################################

      n_ligs = df_extract.shape[0]
      print(n_ligs)
      print(df_extract)
      print(df_extract["canvassim"].unique())

      # group by tranche and extract the list of hits
      df_extract_g = df_extract.groupby('tranche')
      extracted_ligands = []
      for tranche, df_group in df_extract_g:
             print(f'extracting ligands from tranche {tranche}')
             filename = f'{ligpreppdir}//{tranche}_out.maegz'
             ligands = df_group['canvassim'].tolist()
             ligand_indices = df_group['index'].tolist()
             props = {f'r_u_tanimoto_to_{refs[0]}':df_group[refs[0]].to_list()}
              
             print(ligands)
             #extracted_ligands = extracted_ligands + extract_ligs(ligands, filename)
             extracted_ligands = extracted_ligands + extract_ligs_index(ligand_indices, filename,properties=props)


      # save extracted structures
      with structure.StructureWriter(f'{outputfolder}/extracted_ligands.mae') as writer:
               for st in extracted_ligands:
                    print(st.title)
                    writer.append(st)
     # # draw mols
     # rdkit_mols0 = [adapter.to_rdkit(st) for st in extracted_ligands]
     # rdkit_mols = [Chem.RemoveHs(m) for m in rdkit_mols0]
     # for m in rdkit_mols: tmp=AllChem.Compute2DCoords(m)
     # #define the legends
     # propname = f'r_u_tanimoto_to_{refs[0]}'
     # legends=[f'{x.GetProp("_Name")}: \n{float(x.GetProp(propname)):.2f}' for x in rdkit_mols]
     # img=Draw.MolsToGridImage(rdkit_mols,molsPerRow=4,subImgSize=(300,300), legends=legends)
     # img.save('mols.png')  



