
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

from schrodinger import structure
from schrodinger.structutils import analyze


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
      extracted_list = []
      with structure.StructureReader(filename) as reader:
           for st in reader:
                #print(st.title, st.unique_smiles) 
                if st.title in ligand_names:
                       smi = analyze.generate_smiles(st)
                       print(st.title, smi)
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

      #lower_cutoff=0.45
      lower_cutoff=0.4
      ligpreppdir = '/sc/arion/projects/H_filizm02a/ZINC20_3d//20210406_glide_ready_zinc//'

      #read tanimoto data
      tfiles = get_tanimoto_files(outputfolder)
      dff = read_tanimotos(tfiles)
      print(f'all tanimotos read {dff.shape}')

      #extract the names of the reference ligands
      refs = [col for col in dff.columns if not col in ['canvassim','tranche']]
      print(f'reference ligands {refs}')
     
      # sort by the similarity to the first ref ligands
      dff1 =  dff.sort_values(by=refs[0], ascending=False)
      #print(dff1.head())

      df_extract = dff1[dff1[refs[0]] > lower_cutoff]
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
             print(ligands)
             #print("Thanks")
             extracted_ligands = extracted_ligands + extract_ligs(ligands, filename)
      
      tanimoto_scores = []
      for i in extracted_ligands:
             tc_score = df_extract[df_extract["canvassim"]==i.title]["pimethixene"].unique().astype(float)
             #print(tc_score)
             tanimoto_scores.append(tc_score[0])
      print(tanimoto_scores)
             

      # save extracted structures
      with structure.StructureWriter(f'{outputfolder}/extracted_ligands.mae') as writer:
               for st,tc in zip(extracted_ligands,tanimoto_scores):
                    print(st.title)
                    print(tc)
                    st.property["r_user_Tanimoto_Coeff"]=tc
                    writer.append(st)

