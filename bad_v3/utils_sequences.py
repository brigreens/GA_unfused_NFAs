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
import os

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
        e.g. 100_200_230_200_100_01210 for a certain NFA with a specific sequence

    Returns
    -------
    NFA: list (specific format)
        [L_term, core1, core2, core3, R_term, sequence]
    '''

    # splits file name by "_" symbol
    NFA_indices = filename.split('_')
    # creates ordered list of unit index [L_term, core1, core2, core3, R_term, sequence]
    NFA = [NFA_indices[0], NFA_indices[1], NFA_indices[2], NFA_indices[3], NFA_indices[4], str(NFA_indices[5])]

    return NFA

def make_file_name(NFA):
    '''
    Makes file name for a given molecule

    Parameters
    ---------
    NFA: list 
        [L_term_index, core1_index, core2_index, core3_index, R_term_index, sequence]

    Returns
    -------
    file_name: str
        NFA file name (w/o extension) showing unit indicies and sequence
        e.g. 100_200_230_200_100_02120
    '''

    # capture building block indexes as strings for file naming
    term_L = str(NFA[0])
    core1 = str(NFA[1])
    core2 = str(NFA[2])
    core3 = str(NFA[3])
    term_R = str(NFA[4])
    sequence = NFA[5]

    # make file name string for A-D-A'-D-A structure
    file_name = '%s_%s_%s_%s_%s_%s' % (term_L, core1, core2, core3, term_R, sequence)

    return file_name

def make_unit_list():
    '''
    Makes a list of lists of the SMILES of different types of building block units

    Returns
    -------
    unit_list: list
        list of lists of SMILES
        [acceptor_terminal_left, donor_core, acceptor_core, acceptor_terminal_right, donor_terminal_left, donor_terminal_right]
    '''

    don_core = []
    don_term_L = []
    don_term_R = []

    acc_core = []
    acc_term_L = []
    acc_term_R = []

    with open("/ihome/ghutchison/blp62/GA/running_GA/sequence/unit_csvs/acc_core.csv", "r", encoding='utf-8') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        next(csv_reader)
        for row in csv_reader:
            acc_core.append(row[1])

    with open("/ihome/ghutchison/blp62/GA/running_GA/sequence/unit_csvs/acc_left_term_fixed.csv", "r", encoding='utf-8') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        next(csv_reader)
        for row in csv_reader:
            acc_term_L.append(row[1])

    with open("/ihome/ghutchison/blp62/GA/running_GA/sequence/unit_csvs/acc_right_term_fixed.csv", "r", encoding='utf-8') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        next(csv_reader)
        for row in csv_reader:
            acc_term_R.append(row[1])

    with open("/ihome/ghutchison/blp62/GA/running_GA/sequence/unit_csvs/don_core.csv", "r", encoding='utf-8') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        next(csv_reader)
        for row in csv_reader:
            don_core.append(row[1])

    with open("/ihome/ghutchison/blp62/GA/running_GA/sequence/unit_csvs/don_left_term_fixed.csv", "r", encoding='utf-8') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        next(csv_reader)
        for row in csv_reader:
            don_term_L.append(row[1])

    with open("/ihome/ghutchison/blp62/GA/running_GA/sequence/unit_csvs/don_right_term_fixed.csv", "r", encoding='utf-8') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        next(csv_reader)
        for row in csv_reader:
            don_term_R.append(row[1])

    unit_list = [acc_term_L, don_core, acc_core, acc_term_R, don_term_L, don_term_R]

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

def check_orca(filename):
   with open(filename) as f:

        line = f.readline()
        while line:
            if 'aborting the run' in line:
                print(line)
                return False
            line = f.readline()
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

def make_NFA_str(temp_NFA, unit_list):
    '''
    Constructs non-fullerene acceptor (NFA) small molecule SMILES from acceptor and donor building
    blocks SMILES. 

    Parameters
    ---------
    temp_NFA: list
        specific format. temp_NFA = [L_term_index, core1_index, core2_index, core3_index, R_term_index, sequence]

    unit_list: list
        list of lists of strings
        [acceptor_terminal_left, donor_core, acceptor_core, acceptor_terminal_right, donor_terminal_left, donor_terminal_right]

    Returns
    -------
    NFA_string: str
        NFA SMILES string
    '''
    L_termindex = temp_NFA[0]
    core1 = temp_NFA[1]
    core2 = temp_NFA[2]
    core3 = temp_NFA[3]
    R_termindex = temp_NFA[4]
    sequence = list(str(temp_NFA[5]))

    print(sequence)
    # donor left terminal
    if even_or_odd(int(sequence[0])) == 'even':
        NFA_string = unit_list[4][L_termindex]
    # acceptor left terminal
    else:
        NFA_string = unit_list[0][L_termindex]

    # donor core
    if even_or_odd(int(sequence[1])) == 'even':
        NFA_string = NFA_string + unit_list[1][core1]
    # acceptor core
    else:
        NFA_string = NFA_string + unit_list[2][core1]

    # donor core
    if even_or_odd(int(sequence[2])) == 'even':
        NFA_string = NFA_string + unit_list[1][core2]
    # acceptor core
    else:
        NFA_string = NFA_string + unit_list[2][core2]

    # donor core
    if even_or_odd(int(sequence[3])) == 'even':
        NFA_string = NFA_string + unit_list[1][core3]
    # acceptor core
    else:
        NFA_string = NFA_string + unit_list[2][core3]

    # donor right terminal
    if even_or_odd(int(sequence[4])) == 'even':
        NFA_string = NFA_string + unit_list[5][R_termindex]
    # acceptor right terminal
    else:
        NFA_string = NFA_string + unit_list[3][R_termindex]

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
    path_name = '/ihome/ghutchison/blp62/GA/running_GA/sequence/sTDDFT_output/%s.out' % (file_name)

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

def even_or_odd(num):
    '''
    determines if number is even or odd

    Parameters
    --------
    num: int
        number to check if even or odd

    Returns
    -------
    'even' is number is even, 'odd' if number is odd
    '''


    if int(num) % 2 == 0:
        return 'even'
    elif int(num) == 0:
        return 'even'
    else:
        return 'odd'

def update_seq(sequence, position, acc_or_don):
    '''
    Updates 1 position in the sequence, and will fix all numbers following so it goes from lowest to highest #
    '''

    str_seq = list(str(sequence))

    if position == 'core1':
        if acc_or_don == 'acc':
            # current highest acc number
            str_seq[1] = 3
            acc_count = 3

            if even_or_odd(str_seq[2]) == 'odd':
                acc_count +=2
                str_seq[2] = acc_count
                
            if even_or_odd(str_seq[3]) == 'odd':
                acc_count +=2
                str_seq[3] = acc_count

            if even_or_odd(str_seq[4]) == 'odd':
                if str_seq[4] != str_seq[0]:
                    acc_count +=2
                    str_seq[4] = acc_count

        if acc_or_don == 'don':
            # current highest don number
            str_seq[1] = 2
            don_count = 2

            if even_or_odd(str_seq[2]) == 'even':
                don_count +=2
                str_seq[2] = don_count
                
            if even_or_odd(str_seq[3]) == 'even':
                don_count +=2
                str_seq[3] = don_count

            if even_or_odd(str_seq[4]) == 'even':
                if str_seq[4] != str_seq[0]:
                    don_count +=2
                    str_seq[4] = don_count

    if position == 'core2':
        if acc_or_don == 'acc':
            # current highest acc number
            acc_count = 1

            if even_or_odd(str_seq[1]) == 'odd':
                acc_count +=2

            if even_or_odd(str_seq[2]) == 'odd':
                acc_count +=2
                str_seq[2] = acc_count

            if even_or_odd(str_seq[3]) == 'odd':
                acc_count +=2
                str_seq[3] = acc_count

            if even_or_odd(str_seq[4]) == 'odd':
                if str_seq[4] != str_seq[0]:
                    acc_count +=2
                    str_seq[4] = acc_count

        if acc_or_don == 'don':
            don_count = 0
            
            if even_or_odd(str_seq[1]) == 'even':
                don_count +=2

            if even_or_odd(str_seq[2]) == 'even':
                don_count +=2
                str_seq[2] = don_count
                
            if even_or_odd(str_seq[3]) == 'even':
                don_count +=2
                str_seq[3] = don_count

            if even_or_odd(str_seq[4]) == 'even':
                if str_seq[4] != str_seq[0]:
                    don_count +=2
                    str_seq[4] = don_count

    if position == 'core3':
        if acc_or_don == 'acc':
            acc_count = 1

            if even_or_odd(str_seq[1]) == 'odd':
                acc_count +=2

            if even_or_odd(str_seq[2]) == 'odd':
                acc_count +=2

            if even_or_odd(str_seq[3]) == 'odd':
                acc_count +=2
                str_seq[3] = acc_count

            if even_or_odd(str_seq[4]) == 'odd':
                if str_seq[4] != str_seq[0]:
                    acc_count +=2
                    str_seq[4] = acc_count

        if acc_or_don == 'don':
            don_count = 0
            
            if even_or_odd(str_seq[1]) == 'even':
                don_count +=2

            if even_or_odd(str_seq[2]) == 'even':
                don_count +=2
                
            if even_or_odd(str_seq[3]) == 'even':
                don_count +=2
                str_seq[3] = don_count

            if even_or_odd(str_seq[4]) == 'even':
                if str_seq[4] != str_seq[0]:
                    don_count +=2
                    str_seq[4] = don_count

    seq = ''.join(str(e) for e in str_seq)

    return seq

def make_seq_NFA(sequence, smiles_list):
    '''
    creates indices of units to form new NFA

    Parameters
    -----------
    sequence: str
        sequence of units, where odd numbers are acceptors and even numbers are donors
        Ex: '01210' = D-A-D'-A-D structure

    smiles_list: list of lists
        list of all possible unit SMILES with format: [acc_term_L, don_core, acc_core, acc_term_R, don_term_L, don_term_R]

    Returns
    -------
    temp_NFA: list
        specific format: [L_term_index, core1_index, core2_index, core3_index, R_term_index, sequence]
    '''


    str_seq = list(str(sequence))

    position_1 = int(str_seq[0])
    position_2 = int(str_seq[1])
    position_3 = int(str_seq[2])
    position_4 = int(str_seq[3])
    position_5 = int(str_seq[4])

    temp_positions = []
    temp_NFA = []

    # picks the index of smiles for first position (left terminal)
    # donor
    if even_or_odd(position_1) == 'even':
        L_term_index = random.randint(0, len(smiles_list[4])-1)
    # acceptor
    else:
        L_term_index = random.randint(0, len(smiles_list[0])-1)

    temp_NFA.append(L_term_index)
    temp_positions.append(position_1)

    # picks index of smiles for second position (first core)
    if even_or_odd(position_2) == 'even':
        if position_2 in temp_positions:
            temp_index = temp_positions.index(position_2)

            # 923 is the number of donor cores (so as to make sure a terminal is not picked as a core)
            if temp_NFA[temp_index] <= 923:
                temp_NFA.append(temp_NFA[temp_index])
                temp_positions.append('')
            else:
                core1 = random.randint(0, len(smiles_list[1])-1)
                temp_NFA.append(core1)
                temp_positions.append(position_2)
                sequence = update_seq(sequence, 'core1', 'don')
        else:
            core1 = random.randint(0, len(smiles_list[1])-1)
            temp_NFA.append(core1)
            temp_positions.append(position_2)
    else:
        if position_2 in temp_positions:
            temp_index = temp_positions.index(position_2)

            # 233 is the number of acceptor cores (so as to make sure a terminal is not picked as a core)
            if temp_NFA[temp_index] <= 233:
                temp_NFA.append(temp_NFA[temp_index])
                temp_positions.append('')
            else:
                core1 = random.randint(0, len(smiles_list[2])-1)
                temp_NFA.append(core1)
                temp_positions.append(position_2)
                sequence = update_seq(sequence, 'core1', 'acc')
        else:
            core1 = random.randint(0, len(smiles_list[2])-1)
            temp_NFA.append(core1)
            temp_positions.append(position_2)

    # picks index of smiles for third position (second core)
    if even_or_odd(position_3) == 'even':
        if position_3 in temp_positions:
            temp_index = temp_positions.index(position_3)

            # 923 is the number of donor cores (so as to make sure a terminal is not picked as a core)
            if temp_NFA[temp_index] <= 923:
                temp_NFA.append(temp_NFA[temp_index])
                temp_positions.append('')
            else:
                core2 = random.randint(0, len(smiles_list[1])-1)
                temp_NFA.append(core2)
                temp_positions.append(position_3)
                sequence = update_seq(sequence, 'core2', 'don')
        else:
            core2 = random.randint(0, len(smiles_list[1])-1)
            temp_NFA.append(core2)
            temp_positions.append(position_3)
    else:
        if position_3 in temp_positions:
            temp_index = temp_positions.index(position_3)

            # 233 is the number of donor cores (so as to make sure a terminal is not picked as a core)
            if temp_NFA[temp_index] <= 233:
                temp_NFA.append(temp_NFA[temp_index])
                temp_positions.append('')
            else:
                core2 = random.randint(0, len(smiles_list[2])-1)
                temp_NFA.append(core2)
                temp_positions.append(position_3)
                sequence = update_seq(sequence, 'core2', 'acc')
        else:
            core2 = random.randint(0, len(smiles_list[2])-1)
            temp_NFA.append(core2)
            temp_positions.append(position_3)

    # picks index of smiles for fourth position (third core)
    if even_or_odd(position_4) == 'even':
        if position_4 in temp_positions:
            temp_index = temp_positions.index(position_4)

            # 923 is the number of donor cores (so as to make sure a terminal is not picked as a core)
            if temp_NFA[temp_index] <= 923:
                temp_NFA.append(temp_NFA[temp_index])
                temp_positions.append('')
            else:
                core3 = random.randint(0, len(smiles_list[1])-1)
                temp_NFA.append(core3)
                temp_positions.append(position_4)
                sequence = update_seq(sequence, 'core3', 'don')
        else:
            core3 = random.randint(0, len(smiles_list[1])-1)
            temp_NFA.append(core3)
            temp_positions.append(position_4)
    else:
        if position_4 in temp_positions:
            temp_index = temp_positions.index(position_4)

            # 233 is the number of donor cores (so as to make sure a terminal is not picked as a core)
            if temp_NFA[temp_index] <= 233:
                temp_NFA.append(temp_NFA[temp_index])
                temp_positions.append('')
            else:
                core3 = random.randint(0, len(smiles_list[2])-1)
                temp_NFA.append(core3)
                temp_positions.append(position_4)
                sequence = update_seq(sequence, 'core3', 'acc')
        else:
            core3 = random.randint(0, len(smiles_list[2])-1)
            temp_NFA.append(core3)
            temp_positions.append(position_4)

    # picks index of smiles for fifth position (right terminal)
    if even_or_odd(position_5) == 'even':
        if position_5 in temp_positions:
            temp_index = temp_positions.index(position_5)

            # 923 is the number of donor cores (so as to make sure a terminal is not picked as a core)
            if temp_NFA[temp_index] <= 923:
                temp_NFA.append(temp_NFA[temp_index])
                temp_positions.append('')
            else:
                R_term = random.randint(0, len(smiles_list[5])-1)
                temp_NFA.append(R_term)
                temp_positions.append(position_5)
        else:
            R_term = random.randint(0, len(smiles_list[5])-1)
            temp_NFA.append(R_term)
            temp_positions.append(position_5)
    else:
        if position_5 in temp_positions:
            temp_index = temp_positions.index(position_5)

            # 233 is the number of donor cores (so as to make sure a terminal is not picked as a core)
            if temp_NFA[temp_index] <= 233:
                temp_NFA.append(temp_NFA[temp_index])
                temp_positions.append('')
            else:
                R_term = random.randint(0, len(smiles_list[3])-1)
                temp_NFA.append(R_term)
                temp_positions.append(position_5)
        else:
            R_term = random.randint(0, len(smiles_list[3])-1)
            temp_NFA.append(R_term)
            temp_positions.append(position_5)


    temp_NFA.append(sequence)

    return temp_NFA


sym_sequences = ['11111', # A-A-A-A-A
                '11311', # A-A-A'-A-A
                '11011', # A-A-D-A-A
                '13131', # A-A'-A-A'-A
                '13331', # A-A'-A'-A'-A
                '13531', # A-A'-A"-A'-A
                '13031', # A-A'-D-A'-A
                '10101', # A-D-A-D-A
                '10301', # A-D-A'-D-A
                '10001', # A-D-D-D-A
                '10201', # A-D-D'-D-A
                '01010', # D-A-D-A-D
                '01210', # D-A-D'-A-D
                '01110', # D-A-A-A-D
                '01310', # D-A-A'-A-D
                '00000', # D-D-D-D-D
                '00100', # D-D-A-D-D
                '00200', # D-D-D'-D-D
                '02020', # D-D'-D-D'-D
                '02220', # D-D'-D'-D'-D
                '02420', # D-D'-D"-D'-D
                '02120' # D-D'-A-D'-D
        ]


def select_seq():
    '''
    Randomly chooses sequence to arrange acceptor and donor units with

    Returns
    -------
    seq: str
        sequence where even is donor and odd number is acceptor
        Ex: '10211' = A-D-D'-A-A

    '''

    # select symmetric (0) or asymmetric (1)
    sym_or_asym = random.randint(0, 1)

    if sym_or_asym == 0:
        sym_seq_index = random.randint(0, 21)
        seq = sym_sequences[sym_seq_index]


    else: # asymmetric sequence
        index_list = []

        for x in range(5):
            index_list.append(random.randint(0, 9))
        temp_seq = ''.join(str(e) for e in index_list)
        while temp_seq in sym_sequences:
            index_list = []
            for x in range(5):
                index_list.append(random.randint(0, 9))
            temp_seq = ''.join(str(e) for e in index_list)

        # fixes sequence so numbers start from lowest to highest
        new_seq = []
        track_seq = []
        odd_count = 1
        even_count = 0
        
        for x in index_list:
            # checks to see if its a repeat in the sequence
            if x in track_seq:
                temp_index = track_seq.index(x)
                new_seq.append(new_seq[temp_index])
            # if its a new even number, add it and then increase count by 2
            elif even_or_odd(x) == 'even':
                new_seq.append(even_count)
                even_count +=2
            # if its a new odd number, add it and then increase count by 2
            else:
                new_seq.append(odd_count)
                odd_count +=2
            track_seq.append(x)

        seq = ''.join(str(e) for e in new_seq)

        # makes sure not to get stuck in infinite loop
        count_loop = 10
        count = 0
        while seq in sym_sequences:
            # checks to make sure sequence is not symmetrical 
            if count <= count_loop:
                new_seq = []
                track_seq = []
                odd_count = 1
                even_count = 0
                
                for x in index_list:

                    if x in track_seq:
                        temp_index = track_seq.index(x)
                        new_seq.append(new_seq[temp_index])
                    elif even_or_odd(x) == 'even':
                        new_seq.append(even_count)
                        even_count +=2
                    else:
                        new_seq.append(odd_count)
                        odd_count +=2
                    track_seq.append(x)

                seq = ''.join(str(e) for e in new_seq)
                count +=1
            else:
                index_list = []

                for x in range(5):
                    index_list.append(random.randint(0, 9))
                # starts the while loop count over
                count = 0

    return seq

