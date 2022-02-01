import os
import sys
import shutil
import subprocess
import random
from openbabel import pybel
import numpy as np
from scipy import stats
from statistics import mean
from copy import deepcopy
import pickle
from rdkit import Chem
import argparse
import csv
import glob
import math
from itertools import product
import pandas as pd

def make3D(mol):
    '''
    Makes the mol object from SMILES 3D

    Parameters
    ---------
    mol: object
        pybel molecule object
    '''
    # make mol object 3D and add hydrogens
    pybel._builder.Build(mol.OBMol)
    mol.addh()

    ff = pybel._forcefields["mmff94"]
    success = ff.Setup(mol.OBMol)
    if not success:
        ff = pybel._forcefields["uff"]
        success = ff.Setup(mol.OBMol)
        if not success:
            sys.exit("Cannot set up forcefield")

    ff.ConjugateGradients(100, 1.0e-3)
    ff.WeightedRotorSearch(100, 25)
    ff.ConjugateGradients(250, 1.0e-4)

    ff.GetCoordinates(mol.OBMol)


def run_geom_opt(file_name, smiles):
    '''
    Performs geometry optimization with GFN2

    Parameters
    ----------
    NFA: list (specific format)
        [acc_term_index, don_core_index, acc_core_index]
    NFA_str: string
        SMILES string of NFA
    generation_counter: int
        current generation number
    '''

    # make NFA string into pybel molecule object
    mol = pybel.readstring('smi', smiles)
    make3D(mol)

    # write NFA .xyz file to containing folder
    mol.write('xyz', '/ihome/ghutchison/blp62/GA/running_GA/ADDDA/output_files/%s.xyz' % (file_name), overwrite=True)

    # run xTB geometry optimization and sTD-DFT
    calcs = subprocess.call('(cd /ihome/ghutchison/blp62/GA/running_GA/ADDDA/output_files && sbatch -J /ihome/ghutchison/blp62/GA/running_GA/ADDDA/output_files/%s.xyz /ihome/ghutchison/blp62/GA/running_GA/ADDDA/xtb_sTDDFT.slurm)' %(file_name), shell=True)



df = pd.read_csv('/ihome/ghutchison/blp62/GA/running_GA/ADDDA/ADDDA.csv')

smiles = df['smiles']
NFA_names = df['NFA']

for i in range(len(smiles)):
    run_geom_opt(NFA_names[i], smiles[i])