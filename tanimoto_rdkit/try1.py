import os
import warnings
import sys, getopt
import subprocess
import re
import pandas as pd
import numpy as np

from schrodinger import structure
from schrodinger.structutils import analyze
from schrodinger.structutils.analyze import generate_smiles
from schrodinger import adapter

filename = "InducedFit_1-out.maegz"
get_list = []
get_list_name = []
with structure.StructureReader(filename) as reader:
     for i, st in enumerate(reader):
         inf = [i, st.title, st.property['r_i_docking_score']] 
         get_list.append(inf)
         get_list_name.append(st.title)
get_list_name = list(set(get_list_name))
#twolists = [get_list_name, get_list]
#print(get_list)
#print(get_list_name)


get_lowest_energy_list = []
get_lowest_energy_index = []
for i in get_list_name:
    #print(i) 
    lowest_score = 10.0
    for j in get_list:
        if j[1] == i:
               if j[2] < lowest_score:
                    lowest_score = j[2]
                    lowest_score_pose = j
         
    get_lowest_energy_list.append(lowest_score_pose)
    get_lowest_energy_index.append(lowest_score_pose[0])
#twolists =[get_lowest_energy_index, get_lowest_energy_list]
print(get_lowest_energy_list)
print(get_lowest_energy_index)

