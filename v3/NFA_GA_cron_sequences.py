#!/ihome/ghutchison/blp62/.conda/envs/py36/bin/python3


'''
Notes:
-each NFA is defined in the format [acc_term_index, don_core_index, acc_core_index]
-population is list of these NFAs
'''

import os
import sys
import shutil
import subprocess
import random

from numpy.lib.nanfunctions import nanstd
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

import utils_sequences
import scoring_sequences



def run_geom_opt(NFA, NFA_str, gen_counter):
    '''
    Performs geometry optimization with GFN2

    Parameters
    ----------
    NFA: list (specific format)
        [L_term_index, core1_index, core2_index, core3_index, R_term_index, sequence]
    NFA_str: string
        SMILES string of NFA
    generation_counter: int
        current generation number
    '''

    # make file name string w/ convention   
    file_name = utils_sequences.make_file_name(NFA)

    # checks for previous files if it has A-D-A'-D-A structure
    #if NFA[-1] == '10301':
        #orig_filename = file_name.rsplit('_', 1)[0]
        #print(orig_filename)

        # if sTD-DFT exists, add it to generation
    exists = os.path.isfile('/ihome/ghutchison/blp62/GA/running_GA/v5/sTDDFT_output/%s.out' % (file_name))
    if exists:
        duplicate_sTDDFT = subprocess.call('(cp /ihome/ghutchison/blp62/GA/running_GA/v5/sTDDFT_output/%s.out /ihome/ghutchison/blp62/GA/running_GA/v5/generations/%s/%s.out)' % (file_name, gen_counter, file_name), shell=True)
        duplicate_sTDDFT_prop = subprocess.call('(cp /ihome/ghutchison/blp62/GA/running_GA/v5/sTDDFT_output/%s.prop /ihome/ghutchison/blp62/GA/running_GA/v5/generations/%s/%s.prop)' % (file_name, gen_counter, file_name), shell=True)

        print("sTDDFT file existed. Copied to current generation")
        return 
    
    # make NFA string into pybel molecule object
    mol = pybel.readstring('smi', NFA_str)
    utils_sequences.make3D(mol)

    # write NFA .xyz file to containing folder
    mol.write('xyz', '/ihome/ghutchison/blp62/GA/running_GA/v5/generations/%s/%s.xyz' % (gen_counter, file_name), overwrite=True)

    # copy GFN2 input xyz into GFN2_input directory
    duplicate_xtb_inp = subprocess.call('(cp /ihome/ghutchison/blp62/GA/running_GA/v5/generations/%s/%s.xyz /ihome/ghutchison/blp62/GA/running_GA/v5/GFN2_input/%s.xyz)' % (gen_counter, file_name, file_name), shell=True)

    # run xTB geometry optimization and sTD-DFT
    calcs = subprocess.call('(cd /ihome/ghutchison/blp62/GA/running_GA/v5/generations/%s && sbatch -J /ihome/ghutchison/blp62/GA/running_GA/v5/generations/%s/%s.xyz /ihome/ghutchison/blp62/GA/running_GA/v5/xtb_sTDDFT_sequence.slurm)' %(gen_counter, gen_counter, file_name), shell=True)
    

def find_elec_prop(file_names):
    '''
    Calculates oscillator strength and sum of oscillator strength of each NFA in population

    Parameters
    ---------
    file_names: list
        list of NFA filenames (ex: 30_1_9_1_30_10301)

    Returns
    -------
    calc_prop_lists: list
        nested list of [list of difference in HOMO and HOMO-1 energy, list of sum of oscillator strengths]
    '''

    NFA_deltaHOMO_list = []
    NFA_summedoscs_list = []

    # parse sTD-DFT output files
    for NFA in file_names:
        # parses the sTD-DFT output file to find the first oscillator strength and sum of all oscillator strengths 
        stddft_props = utils_sequences.parse_sTDDFT(NFA)
        NFA_deltaHOMO_list.append(stddft_props[0])
        NFA_summedoscs_list.append(stddft_props[1])

    # make nested list of deltaHOMO and summed oscs lists
    stddft_prop_lists = [NFA_deltaHOMO_list, NFA_summedoscs_list]

    return stddft_prop_lists


def parent_select(fitness_list, pop_size):
    '''
    Finds top half of population

    Parameters
    ---------
    fitness_list: list
        fitness_list = [ranked_NFA_names, ranked_PCE, ranked_deltaHOMO, ranked_summed_oscs, ranked_best_donor]
    pop_size: int
        total number of NFAs per generation

    Returns
    -------
    parent_list: list
        top half of population
    '''
    # find number of parents (half of population)
    parent_count = int(pop_size / 2)

    # make list of top NFAs
    parent_list = []
    for x in range(parent_count):
        NFA = []
        for i in range(len(fitness_list[0][x])):
            if i == 5:
                NFA.append(str(fitness_list[0][x][i]))
            else:
                NFA.append(int(fitness_list[0][x][i]))

        parent_list.append(NFA)

    return parent_list


def mutate(NFA, smiles_list, parent_list, parent_b):
    '''
    Creates point mutation in given NFA if NFA is selected for mutation

    Parameters
    ----------
    NFA: list (specific format)
        [L_term_index, core1_index, core2_index, core3_index, R_term_index, sequence]

    smiles_list: nested list (specific format)
        list of monomer SMILES
        [acceptor_terminal_left, donor_core, acceptor_core, acceptor_terminal_right, donor_terminal_left, donor_terminal_right]

    Returns
    -------
    NFA: list (specific format)
        updated NFA with a unit changed
        [L_term_index, core1_index, core2_index, core3_index, R_term_index, sequence]
    '''

    # set mutation rate
    mut_rate = 0.4

    # determine whether to mutate based on mutation rate
    rand = random.randint(1, 10)
    if rand <= (mut_rate * 10):
        pass
    else:
        return NFA

    # choose type of mutation ( 0 = replace unit, 1 = sequence)
    type = random.randint(0, 1)

    if type == 0:

        # get sequence, to make sure replacing units of same kind and preserve sequence
        seq = list(NFA[-1])

        # choose point of mutation (0 = left_term, 1 = core1, 2 = core2, 3 = core3, 4 = right_term)
        point = random.randint(0, 4)

        # replace the left terminal
        if point == 0:
            if utils_sequences.even_or_odd(int(seq[0])) == 'even':
                new_unit = random.randint(0, len(smiles_list[4]) - 1)
            else:
                new_unit = random.randint(0, len(smiles_list[0]) - 1)
            indices = [i for i, x in enumerate(seq) if x == seq[0]]

        # replace core 1
        elif point == 1:
            if utils_sequences.even_or_odd(int(seq[1])) == 'even':
                # donor core
                new_unit = random.randint(0, len(smiles_list[1]) - 1)
            else:
                # acceptor core
                new_unit = random.randint(0, len(smiles_list[2]) - 1)
            indices = [i for i, x in enumerate(seq) if x == seq[1]]

        # replace core 2
        elif point == 2:
            if utils_sequences.even_or_odd(int(seq[2])) == 'even':
                # donor core
                new_unit = random.randint(0, len(smiles_list[1]) - 1)
            else:
                # acceptor core
                new_unit = random.randint(0, len(smiles_list[2]) - 1)
            indices = [i for i, x in enumerate(seq) if x == seq[2]]

        # replace core 3
        elif point == 3:
            if utils_sequences.even_or_odd(int(seq[3])) == 'even':
                # donor core
                new_unit = random.randint(0, len(smiles_list[1]) - 1)
            else:
                # acceptor core
                new_unit = random.randint(0, len(smiles_list[2]) - 1)
            indices = [i for i, x in enumerate(seq) if x == seq[3]]

        # replace right terminal
        elif point == 4:
            if utils_sequences.even_or_odd(int(seq[4])) == 'even':
                # donor core
                new_unit = random.randint(0, len(smiles_list[5]) - 1)
            else:
                # acceptor core
                new_unit = random.randint(0, len(smiles_list[3]) - 1)
            indices = [i for i, x in enumerate(seq) if x == seq[4]]

        for i in indices:
            NFA[i] = new_unit

        return NFA

    # change sequence
    else:
        don_term_units = []
        don_core_units = []
        acc_term_units = []
        acc_core_units = []

        # get sequence, to make sure we can rearrange the unit
        orig_seq_str = str(parent_list[parent_b][-1])
        orig_seq = list(orig_seq_str)

        if utils_sequences.even_or_odd(int(orig_seq[0])) == 'even':
            don_term_units.append(parent_list[parent_b][0])
        else:
            acc_term_units.append(parent_list[parent_b][0])

        if utils_sequences.even_or_odd(int(orig_seq[1])) == 'even':
            if parent_list[parent_b][1] not in don_term_units:
                don_core_units.append(parent_list[parent_b][1])
                don_term_units.append(parent_list[parent_b][1])
        else:
            if parent_list[parent_b][1] not in acc_term_units:
                acc_core_units.append(parent_list[parent_b][1])
                acc_term_units.append(parent_list[parent_b][1])

        if utils_sequences.even_or_odd(int(orig_seq[2])) == 'even':
            if parent_list[parent_b][2] not in don_term_units:
                don_core_units.append(parent_list[parent_b][2])
                don_term_units.append(parent_list[parent_b][2])
        else:
            if parent_list[parent_b][2] not in acc_term_units:
                acc_core_units.append(parent_list[parent_b][2])
                acc_term_units.append(parent_list[parent_b][2])

        if utils_sequences.even_or_odd(int(orig_seq[3])) == 'even':
            if parent_list[parent_b][3] not in don_term_units:
                don_core_units.append(parent_list[parent_b][3])
                don_term_units.append(parent_list[parent_b][3])
        else:
            if parent_list[parent_b][3] not in acc_term_units:
                acc_core_units.append(parent_list[parent_b][3])
                acc_term_units.append(parent_list[parent_b][3])

        if utils_sequences.even_or_odd(int(orig_seq[4])) == 'even':
            if parent_list[parent_b][4] not in don_term_units:
                don_term_units.append(parent_list[parent_b][4])
        else:
            if parent_list[parent_b][4] not in acc_term_units:
                acc_term_units.append(parent_list[parent_b][4])


        # selects new sequence to use
        new_seq_str = str(utils_sequences.select_seq())
        new_seq = list(new_seq_str)
        while new_seq_str == orig_seq_str:
            new_seq_str = str(utils_sequences.select_seq())
            new_seq = list(new_seq_str)

        temp_NFA = []
        
        # track number used so we do not repeat
        don_term_count = 0
        don_core_count = 0
        acc_term_count = 0
        acc_core_count = 0

        seq_type = []

        # new mutation
        # left terminal unit
        if utils_sequences.even_or_odd(int(new_seq[0])) == 'even':
            if len(don_term_units) != 0:
                temp_NFA.append(don_term_units[0])
                don_term_count +=1
            else:
                temp_NFA.append(random.randint(0, len(smiles_list[4])-1))

            seq_type.append('even')
        else:
            if len(acc_term_units) != 0:
                temp_NFA.append(acc_term_units[0])
                acc_term_count +=1
            else:
                temp_NFA.append(random.randint(0, len(smiles_list[0])-1))
            seq_type.append('odd')

        # core 1
        # if same as any previous units
        if int(new_seq[1]) == int(new_seq[0]):
            if utils_sequences.even_or_odd(int(new_seq[0])) == 'even':
                seq_type.append('even')
                # 895 is the number of donor cores (so as to make sure a terminal is not picked as a core)
                if temp_NFA[0] <= 895:
                    temp_NFA.append(temp_NFA[0])
                else:
                    core1 = random.randint(0, len(smiles_list[1])-1)
                    temp_NFA.append(core1)
                    new_seq_str = utils_sequences.update_seq(new_seq_str, 'core1', 'don')
                    new_seq = list(new_seq_str)

            else:
                seq_type.append('odd')
                # 224 is the number of acceptor cores (so as to make sure a terminal is not picked as a core)
                if temp_NFA[0] <= 224:
                    temp_NFA.append(temp_NFA[0])
                else:
                    core1 = random.randint(0, len(smiles_list[2])-1)
                    temp_NFA.append(core1)
                    new_seq_str = utils_sequences.update_seq(new_seq_str, 'core1', 'acc')
                    new_seq = list(new_seq_str)

        else:
            # if new donor
            if utils_sequences.even_or_odd(int(new_seq[1])) == 'even':
                # if there is enough donor units
                if len(don_core_units) > don_core_count:
                    temp_unit = don_core_units[don_core_count]
                    don_core_count +=1
                    # go through past units in sequence to check for duplicates
                    for x in range(len(seq_type)):
                        if seq_type[x] == 'even':
                            while temp_unit == temp_NFA[x]:
                                if len(don_core_units) > don_core_count:
                                    temp_unit = don_core_units[don_core_count]
                                    don_core_count +=1
                                else:
                                    temp_unit = random.randint(0, len(smiles_list[1])-1)
                else:
                    temp_unit = random.randint(0, len(smiles_list[1])-1)
                    # go through past units in sequence to check for duplicates
                    for x in range(len(seq_type)):
                        if seq_type[x] == 'even':
                            while temp_unit == temp_NFA[x]:
                                temp_unit = random.randint(0, len(smiles_list[1])-1)
                temp_NFA.append(temp_unit)
                seq_type.append('even')

            # if new acceptor
            else:
                if len(acc_core_units) > acc_core_count:
                    temp_unit = acc_core_units[acc_core_count]
                    acc_core_count +=1
                    # go through past units in sequence to check for duplicates
                    for x in range(len(seq_type)):
                        if seq_type[x] == 'odd':
                            while temp_unit == temp_NFA[x]:
                                if len(acc_core_units) > acc_core_count:
                                    temp_unit = acc_core_units[acc_core_count]
                                    acc_core_count +=1
                                else:
                                    temp_unit = random.randint(0, len(smiles_list[2])-1)

                else:
                    temp_unit = random.randint(0, len(smiles_list[2])-1)
                    # go through past units in sequence to check for duplicates
                    for x in range(len(seq_type)):
                        if seq_type[x] == 'odd':
                            while temp_unit == temp_NFA[x]:
                                temp_unit = random.randint(0, len(smiles_list[2])-1)
                temp_NFA.append(temp_unit)
                seq_type.append('odd')

        # core 2
        # if same as any previous units
        if int(new_seq[2]) == int(new_seq[0]):
            if utils_sequences.even_or_odd(int(new_seq[0])) == 'even':
                seq_type.append('even')
                # 895 is the number of donor cores (so as to make sure a terminal is not picked as a core)
                if temp_NFA[0] <= 895:
                    temp_NFA.append(temp_NFA[0])
                else:
                    core2 = random.randint(0, len(smiles_list[1])-1)
                    temp_NFA.append(core2)
                    new_seq_str = utils_sequences.update_seq(new_seq_str, 'core2', 'don')
                    new_seq = list(new_seq_str)

            else:
                seq_type.append('odd')
                # 224 is the number of acceptor cores (so as to make sure a terminal is not picked as a core)
                if temp_NFA[0] <= 224:
                    temp_NFA.append(temp_NFA[0])
                else:
                    core2 = random.randint(0, len(smiles_list[2])-1)
                    temp_NFA.append(core2)
                    new_seq_str = utils_sequences.update_seq(new_seq_str, 'core2', 'acc')
                    new_seq = list(new_seq_str)

        elif int(new_seq[2]) == int(new_seq[1]):
            temp_NFA.append(temp_NFA[1])
            if utils_sequences.even_or_odd(int(new_seq[1])) == 'even':
                seq_type.append('even')
            else:
                seq_type.append('odd')
        else:
            # if new donor
            if utils_sequences.even_or_odd(int(new_seq[2])) == 'even':
                # if there is enough donor units
                if len(don_core_units) > don_core_count:
                    temp_unit = don_core_units[don_core_count]
                    don_core_count +=1
                    # go through past units in sequence to check for duplicates
                    for x in range(len(seq_type)):
                        if seq_type[x] == 'even':
                            while temp_unit == temp_NFA[x]:
                                if len(don_core_units) > don_core_count:
                                    temp_unit = don_core_units[don_core_count]
                                    don_core_count +=1
                                else:
                                    temp_unit = random.randint(0, len(smiles_list[1])-1)
                else:
                    temp_unit = random.randint(0, len(smiles_list[1])-1)
                    # go through past units in sequence to check for duplicates
                    for x in range(len(seq_type)):
                        if seq_type[x] == 'even':
                            while temp_unit == temp_NFA[x]:
                                temp_unit = random.randint(0, len(smiles_list[1])-1)
                temp_NFA.append(temp_unit)
                seq_type.append('even')

            # if new acceptor
            else:
                if len(acc_core_units) > acc_core_count:
                    temp_unit = acc_core_units[acc_core_count]
                    acc_core_count +=1
                    # go through past units in sequence to check for duplicates
                    for x in range(len(seq_type)):
                        if seq_type[x] == 'odd':
                            while temp_unit == temp_NFA[x]:
                                if len(acc_core_units) > acc_core_count:
                                    temp_unit = acc_core_units[acc_core_count]
                                    acc_core_count +=1
                                else:
                                    temp_unit = random.randint(0, len(smiles_list[2])-1)
                else:
                    temp_unit = random.randint(0, len(smiles_list[2])-1)
                    # go through past units in sequence to check for duplicates
                    for x in range(len(seq_type)):
                        if seq_type[x] == 'odd':
                            while temp_unit == temp_NFA[x]:
                                temp_unit = random.randint(0, len(smiles_list[2])-1)
                temp_NFA.append(temp_unit)
                seq_type.append('odd')

        # core 3
        # if same as any previous units
        if int(new_seq[3]) == int(new_seq[0]):
            if utils_sequences.even_or_odd(int(new_seq[0])) == 'even':
                seq_type.append('even')
                # 895 is the number of donor cores (so as to make sure a terminal is not picked as a core)
                if temp_NFA[0] <= 895:
                    temp_NFA.append(temp_NFA[0])
                else:
                    core3 = random.randint(0, len(smiles_list[1])-1)
                    temp_NFA.append(core3)
                    new_seq_str = utils_sequences.update_seq(new_seq_str, 'core3', 'don')
                    new_seq = list(new_seq_str)

            else:
                seq_type.append('odd')
                # 224 is the number of acceptor cores (so as to make sure a terminal is not picked as a core)
                if temp_NFA[0] <= 224:
                    temp_NFA.append(temp_NFA[0])
                else:
                    core3 = random.randint(0, len(smiles_list[2])-1)
                    temp_NFA.append(core3)
                    new_seq_str = utils_sequences.update_seq(new_seq_str, 'core3', 'acc')
                    new_seq = list(new_seq_str)

        elif int(new_seq[3]) == int(new_seq[1]):
            temp_NFA.append(temp_NFA[1])
            if utils_sequences.even_or_odd(int(new_seq[1])) == 'even':
                seq_type.append('even')
            else:
                seq_type.append('odd')
        elif int(new_seq[3]) == int(new_seq[2]):
            temp_NFA.append(temp_NFA[2])
            if utils_sequences.even_or_odd(int(new_seq[2])) == 'even':
                seq_type.append('even')
            else:
                seq_type.append('odd')
        else:
            # if new donor
            if utils_sequences.even_or_odd(int(new_seq[3])) == 'even':
                # if there is enough donor units
                if len(don_core_units) > don_core_count:
                    temp_unit = don_core_units[don_core_count]
                    don_core_count +=1
                    # go through past units in sequence to check for duplicates
                    for x in range(len(seq_type)):
                        if seq_type[x] == 'even':
                            while temp_unit == temp_NFA[x]:
                                if len(don_core_units) > don_core_count:
                                    temp_unit = don_core_units[don_core_count]
                                    don_core_count +=1
                                else:
                                    temp_unit = random.randint(0, len(smiles_list[1])-1)
                else:
                    temp_unit = random.randint(0, len(smiles_list[1])-1)
                    # go through past units in sequence to check for duplicates
                    for x in range(len(seq_type)):
                        if seq_type[x] == 'even':
                            while temp_unit == temp_NFA[x]:
                                temp_unit = random.randint(0, len(smiles_list[1])-1)
                temp_NFA.append(temp_unit)
                seq_type.append('even')

            # if new acceptor
            else:
                if len(acc_core_units) > acc_core_count:
                    temp_unit = acc_core_units[acc_core_count]
                    acc_core_count +=1
                    # go through past units in sequence to check for duplicates
                    for x in range(len(seq_type)):
                        if seq_type[x] == 'odd':
                            while temp_unit == temp_NFA[x]:
                                if len(acc_core_units) > acc_core_count:
                                    temp_unit = acc_core_units[acc_core_count]
                                    acc_core_count +=1
                                else:
                                    temp_unit = random.randint(0, len(smiles_list[2])-1)
                else:
                    temp_unit = random.randint(0, len(smiles_list[2])-1)
                    # go through past units in sequence to check for duplicates
                    for x in range(len(seq_type)):
                        if seq_type[x] == 'odd':
                            while temp_unit == temp_NFA[x]:
                                temp_unit = random.randint(0, len(smiles_list[2])-1)
                temp_NFA.append(temp_unit)
                seq_type.append('odd')

        # right terminal
        # if same as any previous units
        if int(new_seq[4]) == int(new_seq[0]):
            temp_NFA.append(temp_NFA[0])
            if utils_sequences.even_or_odd(int(new_seq[0])) == 'even':
                seq_type.append('even')
            else:
                seq_type.append('odd')
        elif int(new_seq[4]) == int(new_seq[1]):
            temp_NFA.append(temp_NFA[1])
            if utils_sequences.even_or_odd(int(new_seq[1])) == 'even':
                seq_type.append('even')
            else:
                seq_type.append('odd')
        elif int(new_seq[4]) == int(new_seq[2]):
            temp_NFA.append(temp_NFA[2])
            if utils_sequences.even_or_odd(int(new_seq[2])) == 'even':
                seq_type.append('even')
            else:
                seq_type.append('odd')
        elif int(new_seq[4]) == int(new_seq[3]):
            temp_NFA.append(temp_NFA[3])
            if utils_sequences.even_or_odd(int(new_seq[3])) == 'even':
                seq_type.append('even')
            else:
                seq_type.append('odd')
        else:
            # if new donor
            if utils_sequences.even_or_odd(int(new_seq[4])) == 'even':
                # if there is enough donor units
                if len(don_term_units) > don_term_count:
                    temp_unit = don_term_units[don_term_count]
                    don_term_count +=1
                    # go through past units in sequence to check for duplicates
                    for x in range(len(seq_type)):
                        if seq_type[x] == 'even':
                            while temp_unit == temp_NFA[x]:
                                if len(don_term_units) > don_term_count:
                                    temp_unit = don_term_units[don_term_count]
                                    don_term_count +=1
                                else:
                                    temp_unit = random.randint(0, len(smiles_list[1])-1)
                else:
                    temp_unit = random.randint(0, len(smiles_list[1])-1)
                    # go through past units in sequence to check for duplicates
                    for x in range(len(seq_type)):
                        if seq_type[x] == 'even':
                            while temp_unit == temp_NFA[x]:
                                temp_unit = random.randint(0, len(smiles_list[1])-1)
                temp_NFA.append(temp_unit)
                seq_type.append('even')

            # if new acceptor
            else:
                if len(acc_term_units) > acc_term_count:
                    temp_unit = acc_term_units[acc_term_count]
                    acc_term_count +=1
                    # go through past units in sequence to check for duplicates
                    for x in range(len(seq_type)):
                        if seq_type[x] == 'odd':
                            while temp_unit == temp_NFA[x]:
                                if len(acc_term_units) > acc_term_count:
                                    temp_unit = acc_term_units[acc_term_count]
                                    acc_term_count +=1
                                else:
                                    temp_unit = random.randint(0, len(smiles_list[2])-1)
                else:
                    temp_unit = random.randint(0, len(smiles_list[2])-1)
                    # go through past units in sequence to check for duplicates
                    for x in range(len(seq_type)):
                        if seq_type[x] == 'odd':
                            while temp_unit == temp_NFA[x]:
                                temp_unit = random.randint(0, len(smiles_list[2])-1)
                temp_NFA.append(temp_unit)
                seq_type.append('odd')


        seq = ''.join(str(e) for e in new_seq)

        temp_NFA.append(seq)

        NFA = temp_NFA

        return NFA


def crossover_mutate(parent_list, smiles_list, current_pop_str):
    '''
    Performs crossover and mutation functions on given population

    Parameters
    ---------
    parent_list: list
        list of parent NFAs
    smiles_list: list
        list of all possible monomer SMILES
        [acceptor_terminal_left, donor_core, acceptor_core, acceptor_terminal_right, donor_terminal_left, donor_terminal_right]
    current_pop_str: list
        running list of SMILES of that generation

    Returns
    -------
    new_pop: list
        population list after crossover and mutation
    '''

    # randomly select two parents (as indexes from parent list) to cross
    parent_a = random.randint(0, len(parent_list) - 1)
    parent_b = random.randint(0, len(parent_list) - 1)

    # ensure parents are unique indiviudals
    if len(parent_list) > 1:
        while parent_b == parent_a:
            parent_b = random.randint(0, len(parent_list) - 1)

    # Obtain sequence from parent_a
    seq = parent_list[parent_a][-1]
    seq_list = list(str(seq))

    don_units = []
    acc_units = []
    for x in seq_list:
        if utils_sequences.even_or_odd(x) == 'even':
            if int(x) not in don_units:
                don_units.append(int(x))
        else:
            if int(x) not in acc_units:
                acc_units.append(int(x))

    don_unit_count = len(don_units)
    acc_unit_count = len(acc_units)

    # obtain types of units from parent b
    temp_seq = parent_list[parent_b][5]
    temp_seq_list = list(str(temp_seq))

    don_term_units = []
    don_core_units = []
    acc_term_units = []
    acc_core_units = []

    if utils_sequences.even_or_odd(int(temp_seq_list[0])) == 'even':
        don_term_units.append(parent_list[parent_b][0])
    else:
        acc_term_units.append(parent_list[parent_b][0])

    if utils_sequences.even_or_odd(int(temp_seq_list[1])) == 'even':
        if parent_list[parent_b][1] not in don_term_units:
            don_core_units.append(parent_list[parent_b][1])
            don_term_units.append(parent_list[parent_b][1])
    else:
        if parent_list[parent_b][1] not in acc_term_units:
            acc_core_units.append(parent_list[parent_b][1])
            acc_term_units.append(parent_list[parent_b][1])

    if utils_sequences.even_or_odd(int(temp_seq_list[2])) == 'even':
        if parent_list[parent_b][2] not in don_term_units:
            don_core_units.append(parent_list[parent_b][2])
            don_term_units.append(parent_list[parent_b][2])
    else:
        if parent_list[parent_b][2] not in acc_term_units:
            acc_core_units.append(parent_list[parent_b][2])
            acc_term_units.append(parent_list[parent_b][2])

    if utils_sequences.even_or_odd(int(temp_seq_list[3])) == 'even':
        if parent_list[parent_b][3] not in don_term_units:
            don_core_units.append(parent_list[parent_b][3])
            don_term_units.append(parent_list[parent_b][3])
    else:
        if parent_list[parent_b][3] not in acc_term_units:
            acc_core_units.append(parent_list[parent_b][3])
            acc_term_units.append(parent_list[parent_b][3])

    if utils_sequences.even_or_odd(int(temp_seq_list[4])) == 'even':
        if parent_list[parent_b][4] not in don_term_units:
            don_term_units.append(parent_list[parent_b][4])
    else:
        if parent_list[parent_b][4] not in acc_term_units:
            acc_term_units.append(parent_list[parent_b][4])


   # create backup lists of units from parent a

    backup_don_term_units = []
    backup_don_core_units = []
    backup_acc_term_units = []
    backup_acc_core_units = []

    if utils_sequences.even_or_odd(int(seq_list[0])) == 'even':
        backup_don_term_units.append(parent_list[parent_a][0])
    else:
        backup_acc_term_units.append(parent_list[parent_a][0])

    if utils_sequences.even_or_odd(int(seq_list[1])) == 'even':
        if parent_list[parent_a][1] not in backup_don_term_units:
            backup_don_core_units.append(parent_list[parent_a][1])
            backup_don_term_units.append(parent_list[parent_a][1])
    else:
        if parent_list[parent_a][1] not in backup_acc_term_units:
            backup_acc_core_units.append(parent_list[parent_a][1])
            backup_acc_term_units.append(parent_list[parent_a][1])

    if utils_sequences.even_or_odd(int(seq_list[2])) == 'even':
        if parent_list[parent_a][2] not in backup_don_term_units:
            backup_don_core_units.append(parent_list[parent_a][2])
            backup_don_term_units.append(parent_list[parent_a][2])
    else:
        if parent_list[parent_a][2] not in backup_acc_term_units:
            backup_acc_core_units.append(parent_list[parent_a][2])
            backup_acc_term_units.append(parent_list[parent_a][2])

    if utils_sequences.even_or_odd(int(seq_list[3])) == 'even':
        if parent_list[parent_a][3] not in backup_don_term_units:
            backup_don_core_units.append(parent_list[parent_a][3])
            backup_don_term_units.append(parent_list[parent_a][3])
    else:
        if parent_list[parent_a][3] not in backup_acc_term_units:
            backup_acc_core_units.append(parent_list[parent_a][3])
            backup_acc_term_units.append(parent_list[parent_a][3])

    if utils_sequences.even_or_odd(int(seq_list[4])) == 'even':
        if parent_list[parent_a][4] not in backup_don_term_units:
            backup_don_term_units.append(parent_list[parent_a][4])
    else:
        if parent_list[parent_a][4] not in backup_acc_term_units:
            backup_acc_term_units.append(parent_list[parent_a][4])

    # create hybrid child
    temp_child = []

    # track number used so we do not repeat
    don_term_count = 0
    backup_don_term_count = 0

    don_core_count = 0
    backup_don_core_count = 0

    acc_term_count = 0
    backup_acc_term_count = 0

    acc_core_count = 0
    backup_acc_core_count = 0

    seq_type = []
    
    # left terminal unit
    if utils_sequences.even_or_odd(int(seq_list[0])) == 'even':
        if len(don_term_units) != 0:
            temp_unit = don_term_units[0]
            don_term_count +=1
        else:
            temp_unit = backup_don_term_units[0]
            backup_don_term_count +=1
        seq_type.append('even')

    else:
        if len(acc_term_units) != 0:
            temp_unit = acc_term_units[0]
            acc_term_count +=1
        else:
            temp_unit = backup_acc_term_units[0]
            backup_acc_term_count +=1
        seq_type.append('odd')

    temp_child.append(temp_unit)


    # core 1
    # if same as any previous units
    if int(seq_list[1]) == int(seq_list[0]):

        if utils_sequences.even_or_odd(int(seq_list[0])) == 'even':
            seq_type.append('even')
            # 895 is the number of donor cores (so as to make sure a terminal is not picked as a core)
            if temp_child[0] <= 895:
                temp_child.append(temp_child[0])
            else:
                core1 = random.randint(0, len(smiles_list[1])-1)
                temp_child.append(core1)
                seq = utils_sequences.update_seq(seq, 'core1', 'don')
                seq_list = list(str(seq))
        else:
            seq_type.append('odd')
            # 224 is the number of acceptor cores (so as to make sure a terminal is not picked as a core)
            if temp_child[0] <= 224:
                temp_child.append(temp_child[0])
            else:
                core1 = random.randint(0, len(smiles_list[2])-1)
                temp_child.append(core1)
                seq = utils_sequences.update_seq(seq, 'core1', 'acc')
                seq_list = list(str(seq))
        
    else:
        # if new donor
        if utils_sequences.even_or_odd(int(seq_list[1])) == 'even':
            # if there is enough donor units
            if len(don_core_units) > don_core_count:
                temp_unit = don_core_units[don_core_count]
                don_core_count +=1
                # go through past units in sequence to check for duplicates
                for x in range(len(seq_type)):
                    if seq_type[x] == 'even':
                        while temp_unit == temp_child[x]:
                            if len(don_core_units) > don_core_count:
                                temp_unit = don_core_units[don_core_count]
                                don_core_count +=1
                            elif len(backup_don_core_units) > backup_don_core_count:
                                temp_unit = backup_don_core_units[backup_don_core_count]
                                backup_don_core_count +=1
                            else:
                                temp_unit = random.randint(0, len(smiles_list[1])-1)

            elif len(backup_don_core_units) > backup_don_core_count:
                temp_unit = backup_don_core_units[don_core_count]
                backup_don_core_count +=1
                # go through past units in sequence to check for duplicates
                for x in range(len(seq_type)):
                    if seq_type[x] == 'even':
                        while temp_unit == temp_child[x]:
                            if len(backup_don_core_units) > backup_don_core_count:
                                temp_unit = backup_don_core_units[backup_don_core_count]
                                backup_don_core_count +=1
                            else:
                                temp_unit = random.randint(0, len(smiles_list[1])-1)
            else:
                temp_unit = random.randint(0, len(smiles_list[1])-1)
            temp_child.append(temp_unit)
            seq_type.append('even')

        # if new acceptor
        else:
            if len(acc_core_units) > acc_core_count:
                temp_unit = acc_core_units[acc_core_count]
                acc_core_count +=1
                # go through past units in sequence to check for duplicates
                for x in range(len(seq_type)):
                    if seq_type[x] == 'odd':
                        while temp_unit == temp_child[x]:
                            if len(acc_core_units) > acc_core_count:
                                temp_unit = acc_core_units[acc_core_count]
                                acc_core_count +=1
                            elif len(backup_acc_core_units) > backup_acc_core_count:
                                temp_unit = backup_acc_core_units[backup_acc_core_count]
                                backup_acc_core_count +=1
                            else:
                                temp_unit = random.randint(0, len(smiles_list[2])-1)
            elif len(backup_acc_core_units) > backup_acc_core_count:
                temp_unit = backup_acc_core_units[acc_core_count]
                backup_acc_core_count +=1
                # go through past units in sequence to check for duplicates
                for x in range(len(seq_type)):
                    if seq_type[x] == 'odd':
                        while temp_unit == temp_child[x]:
                            if len(backup_acc_core_units) > backup_acc_core_count:
                                temp_unit = backup_acc_core_units[backup_acc_core_count]
                                backup_acc_core_count +=1
                            else:
                                temp_unit = random.randint(0, len(smiles_list[2])-1)

            else:
                temp_unit = random.randint(0, len(smiles_list[2])-1)
            temp_child.append(temp_unit)
            seq_type.append('odd')

    # core 2
    # if same as any previous units
    if int(seq_list[2]) == int(seq_list[0]):

        if utils_sequences.even_or_odd(int(seq_list[0])) == 'even':
            seq_type.append('even')
            # 895 is the number of donor cores (so as to make sure a terminal is not picked as a core)
            if temp_child[0] <= 895:
                temp_child.append(temp_child[0])
            else:
                core2 = random.randint(0, len(smiles_list[1])-1)
                temp_child.append(core2)
                seq = utils_sequences.update_seq(seq, 'core2', 'don')
                seq_list = list(str(seq))
        else:
            seq_type.append('odd')
            # 224 is the number of acceptor cores (so as to make sure a terminal is not picked as a core)
            if temp_child[0] <= 224:
                temp_child.append(temp_child[0])
            else:
                core2 = random.randint(0, len(smiles_list[2])-1)
                temp_child.append(core2)
                seq = utils_sequences.update_seq(seq, 'core2', 'acc')
                seq_list = list(str(seq))

    elif int(seq_list[2]) == int(seq_list[1]):
        temp_child.append(temp_child[1])
        if utils_sequences.even_or_odd(int(seq_list[1])) == 'even':
            seq_type.append('even')
        else:
            seq_type.append('odd')
    else:
        # if new donor
        if utils_sequences.even_or_odd(int(seq_list[2])) == 'even':
            # if there is enough donor units
            if len(don_core_units) > don_core_count:
                temp_unit = don_core_units[don_core_count]
                don_core_count +=1
                # go through past units in sequence to check for duplicates
                for x in range(len(seq_type)):
                    if seq_type[x] == 'even':
                        while temp_unit == temp_child[x]:
                            if len(don_core_units) > don_core_count:
                                temp_unit = don_core_units[don_core_count]
                                don_core_count +=1
                            elif len(backup_don_core_units) > backup_don_core_count:
                                temp_unit = backup_don_core_units[backup_don_core_count]
                                backup_don_core_count +=1
                            else:
                                temp_unit = random.randint(0, len(smiles_list[1])-1)
            elif len(backup_don_core_units) > backup_don_core_count:
                temp_unit = backup_don_core_units[don_core_count]
                backup_don_core_count +=1
                # go through past units in sequence to check for duplicates
                for x in range(len(seq_type)):
                    if seq_type[x] == 'even':
                        while temp_unit == temp_child[x]:
                            if len(backup_don_core_units) > backup_don_core_count:
                                temp_unit = backup_don_core_units[backup_don_core_count]
                                backup_don_core_count +=1
                            else:
                                temp_unit = random.randint(0, len(smiles_list[1])-1)
            else:
                temp_unit = random.randint(0, len(smiles_list[1])-1)
            temp_child.append(temp_unit)
            seq_type.append('even')

        # if new acceptor
        else:
            if len(acc_core_units) > acc_core_count:
                temp_unit = acc_core_units[acc_core_count]
                acc_core_count +=1
                # go through past units in sequence to check for duplicates
                for x in range(len(seq_type)):
                    if seq_type[x] == 'odd':
                        while temp_unit == temp_child[x]:
                            if len(acc_core_units) > acc_core_count:
                                temp_unit = acc_core_units[acc_core_count]
                                acc_core_count +=1
                            elif len(backup_acc_core_units) > backup_acc_core_count:
                                temp_unit = backup_acc_core_units[backup_acc_core_count]
                                backup_acc_core_count +=1
                            else:
                                temp_unit = random.randint(0, len(smiles_list[2])-1)
            elif len(backup_acc_core_units) > backup_acc_core_count:
                temp_unit = backup_acc_core_units[backup_acc_core_count]
                backup_acc_core_count +=1
                # go through past units in sequence to check for duplicates
                for x in range(len(seq_type)):
                    if seq_type[x] == 'odd':
                        while temp_unit == temp_child[x]:
                            if len(backup_acc_core_units) > backup_acc_core_count:
                                temp_unit = backup_acc_core_units[backup_acc_core_count]
                                backup_acc_core_count +=1
                            else:
                                temp_unit = random.randint(0, len(smiles_list[2])-1)
            else:
                temp_unit = random.randint(0, len(smiles_list[2])-1)
            temp_child.append(temp_unit)
            seq_type.append('odd')

    # core 3
    # if same as any previous units
    if int(seq_list[3]) == int(seq_list[0]):

        if utils_sequences.even_or_odd(int(seq_list[0])) == 'even':
            seq_type.append('even')
            # 895 is the number of donor cores (so as to make sure a terminal is not picked as a core)
            if temp_child[0] <= 895:
                temp_child.append(temp_child[0])
            else:
                core3 = random.randint(0, len(smiles_list[1])-1)
                temp_child.append(core3)
                seq = utils_sequences.update_seq(seq, 'core3', 'don')
                seq_list = list(str(seq))
        else:
            seq_type.append('odd')
            # 224 is the number of acceptor cores (so as to make sure a terminal is not picked as a core)
            if temp_child[0] <= 224:
                temp_child.append(temp_child[0])
            else:
                core3 = random.randint(0, len(smiles_list[2])-1)
                temp_child.append(core3)
                seq = utils_sequences.update_seq(seq, 'core3', 'acc')
                seq_list = list(str(seq))

    elif int(seq_list[3]) == int(seq_list[1]):
        temp_child.append(temp_child[1])
        if utils_sequences.even_or_odd(int(seq_list[1])) == 'even':
            seq_type.append('even')
        else:
            seq_type.append('odd')
    elif int(seq_list[3]) == int(seq_list[2]):
        temp_child.append(temp_child[2])
        if utils_sequences.even_or_odd(int(seq_list[2])) == 'even':
            seq_type.append('even')
        else:
            seq_type.append('odd')
    else:
        # if new donor
        if utils_sequences.even_or_odd(int(seq_list[3])) == 'even':
            # if there is enough donor units
            if len(don_core_units) > don_core_count:
                temp_unit = don_core_units[don_core_count]
                don_core_count +=1
                # go through past units in sequence to check for duplicates
                for x in range(len(seq_type)):
                    if seq_type[x] == 'even':
                        while temp_unit == temp_child[x]:
                            if len(don_core_units) > don_core_count:
                                temp_unit = don_core_units[don_core_count]
                                don_core_count +=1
                            elif len(backup_don_core_units) > backup_don_core_count:
                                temp_unit = backup_don_core_units[backup_don_core_count]
                                backup_don_core_count +=1
                            else:
                                temp_unit = random.randint(0, len(smiles_list[1])-1)
            elif len(backup_don_core_units) > backup_don_core_count:
                temp_unit = backup_don_core_units[don_core_count]
                backup_don_core_count +=1
                # go through past units in sequence to check for duplicates
                for x in range(len(seq_type)):
                    if seq_type[x] == 'even':
                        while temp_unit == temp_child[x]:
                            if len(backup_don_core_units) > backup_don_core_count:
                                temp_unit = backup_don_core_units[backup_don_core_count]
                                backup_don_core_count +=1
                            else:
                                temp_unit = random.randint(0, len(smiles_list[1])-1)
            else:
                temp_unit = random.randint(0, len(smiles_list[1])-1)
            temp_child.append(temp_unit)
            seq_type.append('even')

        # if new acceptor
        else:
            if len(acc_core_units) > acc_core_count:
                temp_unit = acc_core_units[acc_core_count]
                acc_core_count +=1

                # go through past units in sequence to check for duplicates
                for x in range(len(seq_type)):
                    if seq_type[x] == 'odd':
                        while temp_unit == temp_child[x]:
                            if len(acc_core_units) > acc_core_count:
                                temp_unit = acc_core_units[acc_core_count]
                                acc_core_count +=1
                            elif len(backup_acc_core_units) > backup_acc_core_count:
                                temp_unit = backup_acc_core_units[backup_acc_core_count]
                                backup_acc_core_count +=1
                            else:
                                temp_unit = random.randint(0, len(smiles_list[2])-1)

            elif len(backup_acc_core_units) > backup_acc_core_count:
                temp_unit = backup_acc_core_units[backup_acc_core_count]
                backup_acc_core_count +=1

                # go through past units in sequence to check for duplicates
                for x in range(len(seq_type)):
                    if seq_type[x] == 'odd':

                        while temp_unit == temp_child[x]:
                            if len(backup_acc_core_units) > backup_acc_core_count:

                                temp_unit = backup_acc_core_units[backup_acc_core_count]
                                backup_acc_core_count +=1
                            else:
                                temp_unit = random.randint(0, len(smiles_list[2])-1)
            else:
                temp_unit = random.randint(0, len(smiles_list[2])-1)

            temp_child.append(temp_unit)
            seq_type.append('odd')

    # right terminal
    # if same as any previous units
    if int(seq_list[4]) == int(seq_list[0]):
        temp_child.append(temp_child[0])
        if utils_sequences.even_or_odd(int(seq_list[0])) == 'even':
            seq_type.append('even')
        else:
            seq_type.append('odd')
    elif int(seq_list[4]) == int(seq_list[1]):
        temp_child.append(temp_child[1])
        if utils_sequences.even_or_odd(int(seq_list[1])) == 'even':
            seq_type.append('even')
        else:
            seq_type.append('odd')
    elif int(seq_list[4]) == int(seq_list[2]):
        temp_child.append(temp_child[2])
        if utils_sequences.even_or_odd(int(seq_list[2])) == 'even':
            seq_type.append('even')
        else:
            seq_type.append('odd')
    elif int(seq_list[4]) == int(seq_list[3]):
        temp_child.append(temp_child[3])
        if utils_sequences.even_or_odd(int(seq_list[3])) == 'even':
            seq_type.append('even')
        else:
            seq_type.append('odd')
    else:
        # if new donor
        if utils_sequences.even_or_odd(int(seq_list[4])) == 'even':
            # if there is enough donor units
            if len(don_term_units) > don_term_count:
                temp_unit = don_term_units[don_term_count]
                don_term_count +=1
                # go through past units in sequence to check for duplicates
                for x in range(len(seq_type)):
                    if seq_type[x] == 'even':
                        while temp_unit == temp_child[x]:
                            if len(don_term_units) > don_term_count:
                                temp_unit = don_term_units[don_term_count]
                                don_term_count +=1
                            elif len(backup_don_term_units) > backup_don_term_count:
                                temp_unit = backup_don_term_units[backup_don_term_count]
                                backup_don_term_count +=1
                            else:
                                temp_unit = random.randint(0, len(smiles_list[1])-1)
            elif len(backup_don_term_units) > backup_don_term_count:
                temp_unit = backup_don_term_units[don_term_count]
                backup_don_term_count +=1
                # go through past units in sequence to check for duplicates
                for x in range(len(seq_type)):
                    if seq_type[x] == 'even':
                        while temp_unit == temp_child[x]:
                            if len(backup_don_term_units) > backup_don_term_count:
                                temp_unit = backup_don_term_units[backup_don_term_count]
                                backup_don_term_count +=1
                            else:
                                temp_unit = random.randint(0, len(smiles_list[1])-1)

            else:
                temp_unit = random.randint(0, len(smiles_list[1])-1)
            temp_child.append(temp_unit)
            seq_type.append('even')

        # if new acceptor
        else:
            if len(acc_term_units) > acc_term_count:
                temp_unit = acc_term_units[acc_term_count]
                acc_term_count +=1
                # go through past units in sequence to check for duplicates
                for x in range(len(seq_type)):
                    if seq_type[x] == 'odd':
                        while temp_unit == temp_child[x]:
                            if len(acc_term_units) > acc_term_count:
                                temp_unit = acc_term_units[acc_term_count]
                                acc_term_count +=1
                            elif len(backup_acc_term_units) > backup_acc_term_count:
                                temp_unit = backup_acc_term_units[backup_acc_term_count]
                                backup_acc_term_count +=1
                            else:
                                temp_unit = random.randint(0, len(smiles_list[2])-1)
            elif len(backup_acc_term_units) > backup_acc_term_count:
                temp_unit = backup_acc_term_units[acc_term_count]
                backup_acc_term_count +=1
                # go through past units in sequence to check for duplicates
                for x in range(len(seq_type)):
                    if seq_type[x] == 'odd':
                        while temp_unit == temp_child[x]:
                            if len(backup_acc_term_units) > backup_acc_term_count:
                                temp_unit = backup_acc_term_units[backup_acc_term_count]
                                backup_acc_term_count +=1
                            else:
                                temp_unit = random.randint(0, len(smiles_list[2])-1)

            else: 
                temp_unit = random.randint(0, len(smiles_list[2])-1)
            temp_child.append(temp_unit)
            seq_type.append('odd')
 

    temp_child.append(seq)

    # give child opportunity for mutation
    temp_child = mutate(temp_child, smiles_list, parent_list, parent_b)

    print(temp_child)

    if utils_sequences.even_or_odd(int(seq[1])) == 'even':
        if temp_child[1] > 895:
            return False
    if utils_sequences.even_or_odd(int(seq[1])) == 'odd':
        if temp_child[1] > 224:
            return False
    if utils_sequences.even_or_odd(int(seq[2])) == 'even':
        if temp_child[2] > 895:
            return False
    if utils_sequences.even_or_odd(int(seq[2])) == 'odd':
        if temp_child[2] > 224:
            return False
    if utils_sequences.even_or_odd(int(seq[3])) == 'even':
        if temp_child[3] > 895:
            return False
    if utils_sequences.even_or_odd(int(seq[3])) == 'odd':
        if temp_child[3] > 224:
            return False
    
    temp_child_str = utils_sequences.make_NFA_str(temp_child, smiles_list)

    # prevent duplication in generation
    if temp_child_str in current_pop_str:
        return False


    new_NFA = [temp_child, temp_child_str]
    return new_NFA


def init_gen(pop_size, gen_counter, smiles_list):
    '''
    Initializes parameters, creates population, and runs initial generation

    Parameters
    ----------
    pop_size: int
        number of NFAs in each generation
    gen_counter: int
        current generation number
    smiles_list: list of lists
        list of all possible unit SMILES with format: [acc_term_L, don_core, acc_core, acc_term_R, don_term_L, don_term_R]
    '''

    # set up data directories
    setup = subprocess.call('mkdir /ihome/ghutchison/blp62/GA/running_GA/v5/opt_xyz /ihome/ghutchison/blp62/GA/running_GA/v5/GFN2_input /ihome/ghutchison/blp62/GA/running_GA/v5/GFN2_output /ihome/ghutchison/blp62/GA/running_GA/v5/generations /ihome/ghutchison/blp62/GA/running_GA/v5/sTDDFT_output /ihome/ghutchison/blp62/GA/running_GA/v5/sTDDFT_input', shell=True)

    # create directory for generation 1
    gen1 = subprocess.call('mkdir /ihome/ghutchison/blp62/GA/running_GA/v5/generations/%s' % (gen_counter), shell=True)

    # makes initial unit frequency list
    #utils.init_freq_lists(smiles_list)

    # create inital population as list of NFAs
    population = []
    population_str = []
    counter = 0
    while counter < pop_size:
        # creates list to append the unit indices to
        temp_NFA = []

        # creates sequence of acceptor and donor units. Sequence is a string. Ex: '01220'
        sequence = utils_sequences.select_seq()

        # temp_NFA = [L_term_index, core1_index, core2_index, core3_index, R_term_index, sequence]
        temp_NFA = utils_sequences.make_seq_NFA(sequence, smiles_list)

        # make SMILES string of NFA
        temp_NFA_str = utils_sequences.make_NFA_str(temp_NFA, smiles_list)

        try:
            # convert to canonical SMILES to check for duplication
            canonical_smiles = Chem.MolToSmiles(Chem.MolFromSmiles(temp_NFA_str))

        except:
            # prevents molecules with incorrect valence, which canonical smiles will catch and throw error
            print(temp_NFA_str)
            print('Incorrect valence, could not perform canonical smiles')
            pass

        # check for duplication
        if canonical_smiles in population_str:
            pass
        else:
            population_str.append(canonical_smiles)
            # run xtb and sTD-DFT
            run_geom_opt(temp_NFA, canonical_smiles, gen_counter)

            # updates frequency lists
            #utils.update_frequency_list('acc_term', temp_NFA[0])
            #utils.update_frequency_list('don_core', temp_NFA[1])
            #utils.update_frequency_list('acc_core', temp_NFA[2])
            counter += 1  

    # create new analysis output files
    with open('/ihome/ghutchison/blp62/GA/running_GA/v5/quick_analysis_data.csv', 'w') as quick_file:
        # write analysis file headers
        quick_file.write('gen,min_PCE, max_PCE, med_PCE, min_deltaHOMO,max_deltaHOMO,med_deltaHOMO,min_summedoscs,max_summedoscs,med_summedoscs\n')

    with open('/ihome/ghutchison/blp62/GA/running_GA/v5/full_analysis_data.csv', 'w') as analysis_file:
        # write analysis file headers
        analysis_file.write('gen,filename,deltaHOMO,summedoscs,PCE,donor\n')


def evaluate_prev_gen(gen_counter):
    '''
    Evaluate the output files of the previous generation and determine their fitness

    Parameters
    ----------
    gen_counter: int
        current generation

    Returns
    -------
    fitness_list: list
        list of ranked properties from previous generation
        fitness_list = [ranked_NFA_names, ranked_PCE, ranked_deltaHOMO, ranked_summed_oscs, ranked_best_donor]
    '''

    prev_gen = gen_counter -1
    prev_pop = []
    prev_pop_filename = []


    # removes test.out file
    remove_test = subprocess.call('(cd /ihome/ghutchison/blp62/GA/running_GA/v5/generations/%s && rm /ihome/ghutchison/blp62/GA/running_GA/v5/generations/%s/test.out)' % (prev_gen, prev_gen), shell=True)
    
    for output in glob.iglob('/ihome/ghutchison/blp62/GA/running_GA/v5/generations/%s/*.out' % (prev_gen)):
        
        NFA_filename = output.split('/')[-1].split('.')[0]

        orig_smi = utils_sequences.xyz_to_smiles('/ihome/ghutchison/blp62/GA/running_GA/v5/GFN2_input/%s.xyz' % (NFA_filename))
        # check if geom is good. If not, do not evaluate and pass to next generation
        if utils_sequences.check_geom_opt(orig_smi, '/ihome/ghutchison/blp62/GA/running_GA/v5/opt_xyz/%s.xyz'% (NFA_filename)) == False:
            copy_badgeom = subprocess.call('(cp /ihome/ghutchison/blp62/GA/running_GA/v5/opt_xyz/%s.xyz /ihome/ghutchison/blp62/GA/running_GA/v5/bad_geom_xyz/%s.xyz)' % (NFA_filename, NFA_filename), shell=True)
            
            pass
        elif utils_sequences.check_orca(output) == False:
            pass
        else:

            prev_pop_filename.append(NFA_filename)
        
            # move a copy of the sTD-DFT output files into an overall directory
            #move_files(NFA_filename, prev_gen)

            # turn string of NFA indices into list of indices
            NFA = utils_sequences.make_NFA_from_filename(NFA_filename)
            prev_pop.append(NFA)

    # parses oscillator strength and sum of oscillator strengths from output files
    props_list = find_elec_prop(prev_pop_filename)

    deltaHOMO = props_list[0]
    sum_oscs = props_list[1]

    # find PCE and names of NFA. fitness_list = [ranked_NFA_names, ranked_PCE, ranked_deltaHOMO, ranked_summed_oscs, ranked_best_donor]
    fitness_list = scoring_sequences.PCE_prediction(deltaHOMO, sum_oscs, prev_pop)

    # calculate statistics on the previous generation
    min_PCE = fitness_list[1][-1]
    max_PCE = fitness_list[1][0]
    median = int((len(fitness_list[1])-1)/2)
    med_PCE = fitness_list[1][median]

    min_test_deltaHOMO = fitness_list[2][-1][0]
    max_test_deltaHOMO = fitness_list[2][0][0]
    med_test_deltaHOMO = fitness_list[2][median][0]

    min_test_summedoscs = fitness_list[3][-1]
    max_test_summedoscs = fitness_list[3][0]
    med_test_summedoscs = fitness_list[3][median]
 
    with open('/ihome/ghutchison/blp62/GA/running_GA/v5/quick_analysis_data.csv', 'a') as quick_file:
        # write to quick analysis file
        quick_file.write('%d,%f,%f,%f,%f,%f,%f,%f,%f,%f\n' % (prev_gen, min_PCE, max_PCE, med_PCE, min_test_deltaHOMO, max_test_deltaHOMO, med_test_deltaHOMO, min_test_summedoscs, max_test_summedoscs, med_test_summedoscs))

    # find values and write to full analysis file
    # loop over every NFA in population
    for x in range(len(fitness_list[0])):
        NFA = fitness_list[0][x]
        file_name = utils_sequences.make_file_name(NFA)
        donor = fitness_list[4][x]
        PCE = fitness_list[1][x]
        deltaHOMO = fitness_list[2][x][0]
        summed_oscs = fitness_list[3][x]
        
        with open('/ihome/ghutchison/blp62/GA/running_GA/v5/full_analysis_data.csv', 'a') as analysis_file:
            analysis_file.write('%d,%s,%f,%f,%f,%s,\n' % (prev_gen, file_name, deltaHOMO, summed_oscs, PCE, donor))

    # make backup copies of output files
    shutil.copy('/ihome/ghutchison/blp62/GA/running_GA/v5/quick_analysis_data.csv', '/ihome/ghutchison/blp62/GA/running_GA/v5/quick_analysis_data_copy.csv')
    shutil.copy('/ihome/ghutchison/blp62/GA/running_GA/v5/full_analysis_data.csv', '/ihome/ghutchison/blp62/GA/running_GA/v5/full_analysis_data_copy.csv')

    # fitness_list = [ranked_NFA_names, ranked_PCE, ranked_deltaHOMO, ranked_summed_oscs, ranked_best_donor]
    return fitness_list


def next_gen(gen_counter, pop_size, smiles_list, fitness_list):
    '''
    Runs the next generation

    Parameters
    ---------
    gen_counter: int
        current generation running
    pop_size: int
        number of NFAs per generation
    smiles_list: list
        nested list of SMILES
        [acceptor_terminal_left, donor_core, acceptor_core, acceptor_terminal_right, donor_terminal_left, donor_terminal_right]
    fitness_list: list
        fitness_list = [ranked_NFA_names, ranked_PCE, ranked_deltaHOMO, ranked_summed_oscs, ranked_best_donor]
    '''

    # makes a directory in generations directory with current genetation number
    setup = subprocess.call('mkdir /ihome/ghutchison/blp62/GA/running_GA/v5/generations/%s' % (gen_counter), shell=True)

    # Selection - select heaviest (best) 50% of NFAs as parents
    parents = parent_select(fitness_list, pop_size)

    new_pop = deepcopy(parents)
    new_pop_str = []
    for parent in new_pop:
        parent_str = utils_sequences.make_NFA_str(parent, smiles_list)
        new_pop_str.append(parent_str)
        # make file name string w/ convention   
        file_name = utils_sequences.make_file_name(parent)
        duplicate_sTDDFT = subprocess.call('(cp /ihome/ghutchison/blp62/GA/running_GA/v5/sTDDFT_output/%s.out /ihome/ghutchison/blp62/GA/running_GA/v5/generations/%s/%s.out)' % (file_name, gen_counter, file_name), shell=True)
        duplicate_sTDDFT_prop = subprocess.call('(cp /ihome/ghutchison/blp62/GA/running_GA/v5/sTDDFT_output/%s.prop /ihome/ghutchison/blp62/GA/running_GA/v5/generations/%s/%s.prop)' % (file_name, gen_counter, file_name), shell=True)
        #run_sTDDFT(gen_counter, parent)


    while len(new_pop) < pop_size:

        cross_mut_pop = crossover_mutate(parents, smiles_list, new_pop_str)
        # checks for duplication in generation
        while cross_mut_pop == False:
            cross_mut_pop = crossover_mutate(parents, smiles_list, new_pop_str)
        
        new_NFA = cross_mut_pop[0]
        new_NFA_str = cross_mut_pop[1]

        run_geom_opt(new_NFA, new_NFA_str, gen_counter)

        new_pop.append(new_NFA)
        new_pop_str.append(new_NFA_str)

        #filename = utils.make_file_name(new_NFA)

        #utils.update_frequency_list('acc_term', new_NFA[0])
        #utils.update_frequency_list('don_core', new_NFA[1])
        #utils.update_frequency_list('acc_core', new_NFA[2])



def main(gen_counter):
    '''
    main program to run the GA

    Parameters
    ---------
    gen_counter: int
        generation to run next
    '''

    # number of NFAs in population
    pop_size = 32

    # Create list of possible building block unit SMILES in specific format ([acc_term_L, don_core, acc_core, acc_term_R, don_term_L, don_term_R])
    smiles_list = utils_sequences.make_unit_list()
            
    # run if this is the initial generation
    if gen_counter == 1:
        # initialize population and submit calculations
        init_gen(pop_size, gen_counter, smiles_list)
        
        # pickle random state for restart
        randstate = random.getstate()
        rand_file = open('/ihome/ghutchison/blp62/GA/running_GA/v5/rand_states/randstate_%s.p' % (str(gen_counter)), 'wb')
        pickle.dump(randstate, rand_file)
        rand_file.close()

    # run if NOT the initial generation
    else:
        # reload random state

        prev_gen = str(gen_counter -1)

        open_rand = open('/ihome/ghutchison/blp62/GA/running_GA/v5/rand_states/randstate_%s.p' % (prev_gen), 'rb')
        randstate = pickle.load(open_rand)
        random.setstate(randstate)
        open_rand.close()

        # evaluate previous generation
        fitness_list = evaluate_prev_gen(gen_counter)

        # run next generation of GA
        next_gen(gen_counter, pop_size, smiles_list, fitness_list)

        # save pickle random state for restart
        randstate = random.getstate()
        rand_file = open('/ihome/ghutchison/blp62/GA/running_GA/v5/rand_states/randstate_%s.p' % (str(gen_counter)), 'wb')
        pickle.dump(randstate, rand_file)
        rand_file.close()

          
def check_prev_gen_complete(gen_counter):
    '''
    checks to see if the previous generation is complete. If so, it will run the next generation

    Parameters
    ----------
    gen_counter (int)
        generation to run nxt
    '''

    # TODO: add to check if all geoms are good. If they are not good, need to count how many molecules and replace it

    if gen_counter == 1:
        return True
    else:
        prev_gen = gen_counter -1
        pop_size = 32

        # checks to see if all of the sTDDFT calculations are finished. *.prop used since it only appears when calculation is finished. 
        prop_counter = 0
        for x in glob.iglob('/ihome/ghutchison/blp62/GA/running_GA/v5/generations/%s/*.prop' % (prev_gen)):
            prop_counter +=1
        if prop_counter == pop_size:
            return True
        else:
            return False
            

if __name__ == '__main__':
    usage = "usage: %prog [options] "
    parser = argparse.ArgumentParser(usage)

    # sets gen input argument
    parser.add_argument('--gen', action='store', type=int, required = True)
    args = parser.parse_args()

    # total number of generations you would like to run
    max_gen = args.gen
    
    # count number of subdirectories currently in generations directory 
    rootdir = '/ihome/ghutchison/blp62/GA/running_GA/v5/generations/'
    gen_counter = 0
    for dir in os.listdir(rootdir):
        full_path = rootdir + dir
        # count only directories (not files) and ignore hidden directories
        if (os.path.isdir(full_path)):
            gen_counter += 1
    # add one so gen_counter shows next generation to be run
    gen_counter += 1
    print(gen_counter)
    
    # checks to make sure all of the previous generation's sTD-DFT files are complete
    if check_prev_gen_complete(gen_counter) != True:
        sys.exit()
    # runs next generation
    elif gen_counter <= max_gen:      
        main(gen_counter)
    # quits if all generations have been run
    else:
        sys.exit()
        
        
        