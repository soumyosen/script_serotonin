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



