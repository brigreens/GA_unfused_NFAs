'''
Utility script for general functions such as making
smiles strings and files
'''
from openbabel import pybel
import random
import math
from copy import deepcopy
from itertools import product
import csv
import pandas as pd
import subprocess
import numpy as np

ob = pybel.ob


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

def make_NFA_from_filename(filename):
    '''
    takes filename and converts it back to NFA list structure
    
    Parameters
    ----------
    filename: str
        NFA file name (w/o extension) showing unit indicies and full sequence
        e.g. 100_200_230_200_100 for a certain NFA with A-D-A'-D-A structure

    Returns
    -------
    NFA: list (specific format)
        [acc_term_index, don_core_index, acc_core_index]
    '''

    # splits file name by "_" symbol
    NFA_indices = filename.split('_')
    # creates ordered list of unit index [acc_term, don_core, acc_core]
    NFA = [NFA_indices[0], NFA_indices[1], NFA_indices[2]]

    return NFA

def make_file_name(NFA):
    '''
    Makes file name for a given molecule

    Parameters
    ---------
    NFA: list 
        [terminal acceptor index, core acceptor index, core donor index]

    Returns
    -------
    file_name: str
        NFA file name (w/o extension) showing unit indicies
        e.g. 100_200_230_200_100 for a certain NFA with A-D-A'-D-A structure
    '''

    # capture building block indexes as strings for file naming
    termA = str(NFA[0])
    coreD = str(NFA[1])
    coreA = str(NFA[2])

    # make file name string for A-D-A'-D-A structure
    file_name = '%s_%s_%s_%s_%s' % (termA, coreD, coreA, coreD, termA)

    return file_name

def make_unit_list():
    '''
    Makes a list of lists of the SMILES of different types of building block units

    Returns
    -------
    unit_list: list
        list of lists of SMILES
        [acceptor_terminal_left, donor_core, acceptor_core, acceptor_terminal_right]
    '''

    acc_term_L = []
    don_core = []
    acc_core = []
    acc_term_R = []

    with open("/ihome/ghutchison/blp62/GA/running_GA/more_units/unit_csvs/acc_core.csv", "r") as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        next(csv_reader)
        for row in csv_reader:
            acc_core.append(row[1])

    with open("/ihome/ghutchison/blp62/GA/running_GA/more_units/unit_csvs/acc_left_term.csv", "r") as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        next(csv_reader)
        for row in csv_reader:
            acc_term_L.append(row[1])

    with open("/ihome/ghutchison/blp62/GA/running_GA/more_units/unit_csvs/acc_right_term.csv", "r") as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        next(csv_reader)
        for row in csv_reader:
            acc_term_R.append(row[1])

    with open("/ihome/ghutchison/blp62/GA/running_GA/more_units/unit_csvs/don_core.csv", "r") as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        next(csv_reader)
        for row in csv_reader:
            don_core.append(row[1])

    unit_list = [acc_term_L, don_core, acc_core, acc_term_R]

    return unit_list

def xyz_to_smiles(filename):
    '''
    converts an xyz file into SMILES string

    Parameters
    ----------
    filename: str
        file path to xyz file

    Returns
    -------
    SMILES string of that molecule
    '''
    
    mol = next(pybel.readfile("xyz", filename))

    smi = mol.write(format="smi")

    return smi.split()[0].strip()

def RepresentsInt(s):
    '''
    Checks if character is an integer

    Parameters
    ----------
    s: str
        character in string

    Returns
    -------
    True if it is an integer, False if it is not
    '''
    try: 
        int(s)
        return True
    except ValueError:
        return False

def num_rings(SMILES):
    '''
    Counts the number of rings in SMILES

    Parameters
    ----------
    SMILES: str
        SMILES string of the molecule

    Returns
    --------
    ring_count: int
        number of rings in the molecule

    '''

    ring_symbol_count = 0
    frags = list(SMILES)

    for x in frags:
        if RepresentsInt(x) == True:
            ring_symbol_count += 1 # adds 1 if it sees a number in the SMILES, representing a part of the ring break
        elif x == '%':
            ring_symbol_count -= 1 # the number of the ring closure is double digits, so we don't want to add that ring twice

    ring_count = ring_symbol_count / 2 # needs 2 numbers for every 1 ring break

    return ring_count

def check_mol_breaks(SMILES):
    '''
    Checks to see if the molecule was broken into fragments during geometry optimization

    Parameters
    ----------
    SMILES:str
        SMILES string of the molecule

    '''
    frags = list(SMILES)
    if '.' in frags:
        return True

def check_geom_opt(NFA_str, file_name):
    '''
    Checks to see if something weird and incorrect happened during geometry optimization
    Primary example is new rings forming or fragments breaking

    Parameters
    ----------
    NFA_str: string
        SMILES string of the molecule

    file_name: string
        path to xyz file

    '''

    unopt_num_rings = num_rings(NFA_str)

    opt_smi = xyz_to_smiles(file_name)
    opt_num_rings = num_rings(opt_smi)

    if unopt_num_rings != opt_num_rings:
        print(file_name)
        print('The numbers of rings does match before and after geometry optimization')
        return False

    if check_mol_breaks(opt_smi) == True:
        print(file_name)
        print('The molecule broke into fragments')
        return False

def make_NFA_str(acc_term, don_core, acc_core, unit_list):
    '''
    Constructs non-fullerene acceptor (NFA) small molecule SMILES from acceptor and donor building
    blocks SMILES. Structure is A-D-A'-D-A

    Parameters
    ---------
    acc_term: int
        index of acceptor terminal unit
    don_core: int
        index of donor core unit
    acc_core: int
        index of acceptor core unit
    unit_list: list
        list of lists of strings
        [acceptor_terminal_left, donor_core, acceptor_core, acceptor_terminal_right]

    Returns
    -------
    NFA_string: str
        NFA SMILES string
    '''

    NFA_string = unit_list[0][acc_term] + unit_list[1][don_core] + unit_list[2][acc_core]+ unit_list[1][don_core]+ unit_list[3][acc_term]

    return NFA_string


def parse_sTDDFT(file_name):
    '''
    Parses through sTD-DFT output files

    Parameters
    -----------
    file_name: str
        path to output file

    Returns
    -------
    stddft_prop: list
        returns list of difference in energy between HOMO and HOMO-1 and sum of all oscillator strengths
    '''
    path_name = '/ihome/ghutchison/blp62/GA/running_GA/sTDDFT_output/%s.out' % (file_name)

    with open(path_name, 'r', encoding = 'utf-8') as file:
        line = file.readline()
        oscs = []
        deltaHOMO = []
        while line:
            if 'ordered frontier orbitals' in line:
                for x in range(11):
                    line = file.readline()
                HOMOminus1 = float(line[9:15])
                line = file.readline()
                HOMO = float(line[9:15])
                deltaHOMO.append(abs(HOMOminus1 - HOMO))

            if 'ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS' in line:
                for x in range(5):
                    line = file.readline()                
                while line != '\n':
                    oscs.append(float(line[25:37]))
                    line = file.readline()
            line = file.readline()  
        line = file.readline()
   
    summed_oscs = np.sum(oscs)

    stddft_prop = [deltaHOMO, summed_oscs]

    return stddft_prop

def init_freq_lists(smiles_list):
    '''
    Creates csv files with all of the units, and sets the frequency (number of times it appears in an NFA) to 0

    Parameters
    ---------
    smiles_list: list (specific format)
        list of lists of strings
        [acceptor_terminal_left, donor_core, acceptor_core, acceptor_terminal_right]

    '''
    #freq_dir = subprocess.call('mkdir frequency_lists', shell=True)

    with open('/ihome/ghutchison/blp62/GA/running_GA/more_units/frequency_lists/freq_acc_term.csv', 'w') as freq_acc_term:
        fieldnames = ['index', 'SMILES', 'frequency']
        writer = csv.DictWriter(freq_acc_term, fieldnames = fieldnames)
        writer.writeheader()
        
        for acc_term_index in range(len(smiles_list[0])):
            writer.writerow({'index': acc_term_index, 'SMILES': smiles_list[0][acc_term_index], 'frequency': 0})

    with open('/ihome/ghutchison/blp62/GA/running_GA/more_units/frequency_lists/freq_acc_core.csv', 'w') as freq_acc_core:
        fieldnames = ['index', 'SMILES', 'frequency']
        writer = csv.DictWriter(freq_acc_core, fieldnames = fieldnames)
        writer.writeheader()
        
        for acc_core_index in range(len(smiles_list[2])):
            writer.writerow({'index': acc_core_index,'SMILES': smiles_list[2][acc_core_index], 'frequency': 0})

    with open('/ihome/ghutchison/blp62/GA/running_GA/more_units/frequency_lists/freq_don_core.csv', 'w') as freq_don_core:
        fieldnames = ['index', 'SMILES', 'frequency']
        writer = csv.DictWriter(freq_don_core, fieldnames = fieldnames)
        writer.writeheader()
        
        for don_core_index in range(len(smiles_list[1])):
            writer.writerow({'index': don_core_index, 'SMILES': smiles_list[1][don_core_index],'frequency': 0})


def update_frequency_list(unit_type, index):
    '''
    updates the frequency csv files if the unit is used again in an NFA

    Parameters
    ---------
    unit_type: str
        options are 'acc_term', 'acc_core', or 'don_core'
    index: int
        index of the building block unit
    '''

    if unit_type == 'acc_term':
        unit_csv = "/ihome/ghutchison/blp62/GA/running_GA/more_units/frequency_lists/freq_acc_term.csv"

    elif unit_type =='acc_core':
        unit_csv = "/ihome/ghutchison/blp62/GA/running_GA/more_units/frequency_lists/freq_acc_core.csv"

    elif unit_type == 'don_core':
        unit_csv = "/ihome/ghutchison/blp62/GA/running_GA/more_units/frequency_lists/freq_don_core.csv"
    
    # reading the csv
    df = pd.read_csv(unit_csv)

    # updating frequency of unit index
    df.loc[index, 'frequency'] += 1

    # writing into the file
    df.to_csv(unit_csv, index=False)

