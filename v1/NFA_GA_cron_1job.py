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

import utils
import scoring

def run_geom_opt(NFA, NFA_str, gen_counter):
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

    # make file name string w/ convention   
    file_name = utils.make_file_name(NFA)

    # if sTD-DFT exists, add it to generation
    exists = os.path.isfile('/ihome/ghutchison/blp62/GA/running_GA/sTDDFT_output/%s.out' % (file_name))
    if exists:
        duplicate_sTDDFT = subprocess.call('(cp /ihome/ghutchison/blp62/GA/running_GA/sTDDFT_output/%s.out /ihome/ghutchison/blp62/GA/running_GA/generations/%s/%s.out)' % (file_name, gen_counter, file_name), shell=True)
        duplicate_sTDDFT_prop = subprocess.call('(cp /ihome/ghutchison/blp62/GA/running_GA/sTDDFT_output/%s.prop /ihome/ghutchison/blp62/GA/running_GA/generations/%s/%s.prop)' % (file_name, gen_counter, file_name), shell=True)
        print("sTDDFT file existed. Copied to current generation")
        return 

    else:
        # make NFA string into pybel molecule object
        mol = pybel.readstring('smi', NFA_str)
        utils.make3D(mol)

        # write NFA .xyz file to containing folder
        mol.write('xyz', '/ihome/ghutchison/blp62/GA/running_GA/generations/%s/%s.xyz' % (gen_counter, file_name), overwrite=True)

        # copy GFN2 input xyz into GFN2_input directory
        duplicate_xtb_inp = subprocess.call('(cp /ihome/ghutchison/blp62/GA/running_GA/generations/%s/%s.xyz /ihome/ghutchison/blp62/GA/running_GA/GFN2_input/%s.xyz)' % (gen_counter, file_name, file_name), shell=True)

        # run xTB geometry optimization and sTD-DFT
        calcs = subprocess.call('(cd /ihome/ghutchison/blp62/GA/running_GA/generations/%s && sbatch -J /ihome/ghutchison/blp62/GA/running_GA/generations/%s/%s.xyz /ihome/ghutchison/blp62/xtb_sTDDFT.slurm)' %(gen_counter, gen_counter, file_name), shell=True)


def run_sTDDFT(gen_counter, NFA):
    '''
    Submits NFAs for sTD-DFT calculation with slurm
    ---------
    gen_counter: int
        current generation
    NFA: list (specific format)
        [acc_term_index, don_core_index, acc_core_index]
    '''

    # make file name string w/ convention   
    file_name = utils.make_file_name(NFA)

    # if sTD-DFT exists, add it to generation
    exists = os.path.isfile('/ihome/ghutchison/blp62/GA/running_GA/sTDDFT_output/%s.out' % (file_name))
    if exists:
        duplicate_sTDDFT = subprocess.call('(cp /ihome/ghutchison/blp62/GA/running_GA/sTDDFT_output/%s.out /ihome/ghutchison/blp62/GA/running_GA/generations/%s/%s.out)' % (file_name, gen_counter, file_name), shell=True)
        duplicate_sTDDFT_prop = subprocess.call('(cp /ihome/ghutchison/blp62/GA/running_GA/sTDDFT_output/%s.prop /ihome/ghutchison/blp62/GA/running_GA/generations/%s/%s.prop)' % (file_name, gen_counter, file_name), shell=True)

        print("sTDDFT file existed. Copied to current generation")
        return 
    else:
        # move optimized xyz file into sTDDFT input folder
        create_input = subprocess.call('(cd /ihome/ghutchison/blp62/GA/running_GA/generations/%s && obabel /ihome/ghutchison/blp62/GA/running_GA/opt_xyz/%s.xyz -O .inp -o orcainp -m -xf /ihome/ghutchison/blp62/GA/running_GA/cam_b3lyp_def2_SVP_pp.txt)' % (gen_counter, file_name), shell=True)
        # submit sTD-DFT calculation
        sTDDFT_submit = subprocess.call('(cd /ihome/ghutchison/blp62/GA/running_GA/generations/%s && sbatch -J %s.inp /ihome/ghutchison/blp62/orca4_sTDDFT_pp.slurm)' % (gen_counter, file_name), shell=True)

def move_files(file_name, prev_gen):
    '''
    moves output files from previous generation into organized directories

    Parameters
    ----------
    file_name: str
        NFA file name (w/o extension) showing unit indicies
        e.g. 100_200_230_200_100 for a certain NFA with A-D-A'-D-A structure
    prev_gen: int
        Number of the previous generation
    '''

    # makes a copy of the sTDDFT output and input files into separate directories
    copy_sTDDFT_output = subprocess.call('(cp /ihome/ghutchison/blp62/GA/running_GA/generations/%s/%s.out /ihome/ghutchison/blp62/GA/running_GA/sTDDFT_output/%s.out)' % (prev_gen, file_name, file_name), shell=True)
    copy_sTDDFT_prop = subprocess.call('(cp /ihome/ghutchison/blp62/GA/running_GA/generations/%s/%s.prop /ihome/ghutchison/blp62/GA/running_GA/sTDDFT_output/%s.prop)' % (prev_gen, file_name, file_name), shell=True)


def find_elec_prop(file_names):
    '''
    Calculates oscillator strength and sum of oscillator strength of each NFA in population

    Parameters
    ---------
    file_names: list
        list of NFA filenames (ex: 30_1_9_1_30)

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
        stddft_props = utils.parse_sTDDFT(NFA)
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
            NFA.append(int(fitness_list[0][x][i]))

        parent_list.append(NFA)

    return parent_list


def mutate(NFA, smiles_list):
    '''
    Creates point mutation in given NFA if NFA is selected for mutation

    Parameters
    ----------
    NFA: list (specific format)
        [acc_term_index, don_core_index, acc_core_index]

    smiles_list: nested list (specific format)
        list of monomer SMILES
        [acc_term_left_list, don_core_list, acc_core_list, acc_term_right_list]

    Returns
    -------
    NFA: list (specific format)
        updated NFA with a unit changed
        [acc_term_index, don_core_index, acc_core_index]
    '''

    # set mutation rate
    mut_rate = 0.4

    # determine whether to mutate based on mutation rate
    rand = random.randint(1, 10)
    if rand <= (mut_rate * 10):
        pass
    else:
        return NFA

    # choose point of mutation (0 = acc_term, 1 = don_core, 2 = acc_core)
    point = random.randint(0, 2)

    # replace the acceptor terminal
    if point == 0:
        new_unit = random.randint(0, len(smiles_list[0]) - 1)

    # replace the donor core
    elif point == 1:
        new_unit = random.randint(0, len(smiles_list[1]) - 1)

    # replace the acceptor core
    else:
        new_unit = random.randint(0, len(smiles_list[2]) - 1)

    NFA[point] = new_unit

    return NFA


def crossover_mutate(parent_list, pop_size, smiles_list, current_pop_str):
    '''
    Performs crossover and mutation functions on given population

    Parameters
    ---------
    parent_list: list
        list of parent NFAs
    pop_size: int
        number of NFAs in each generation
    smiles_list: list
        list of all possible monomer SMILES
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

    # determine which units to swap. 
    type_unit = random.randint(0, 2)

    # create hybrid child
    temp_child = []

    # 0: acc_term from parent a, don_core and acc core from parent b
    if type_unit == 0:
        temp_child.append(parent_list[parent_a][0])
        temp_child.append(parent_list[parent_b][1])
        temp_child.append(parent_list[parent_b][2])

    # 1: don_core from parent a, acc_term and acc core from parent b
    elif type_unit == 1:
        temp_child.append(parent_list[parent_b][0])
        temp_child.append(parent_list[parent_a][1])
        temp_child.append(parent_list[parent_b][2])

    # 2: acc_core from parent a, don_core and acc_term from parent b
    elif type_unit == 2:
        temp_child.append(parent_list[parent_b][0])
        temp_child.append(parent_list[parent_b][1])
        temp_child.append(parent_list[parent_a][2])

    # give child opportunity for mutation
    temp_child = mutate(temp_child, smiles_list)

    temp_child_str = utils.make_NFA_str(temp_child[0], temp_child[1], temp_child[2], smiles_list)

    # prevent duplication in generation
    if temp_child_str in current_pop_str:
        return False
    #else:
        # updates the unit frequency list for the new NFAs
        #utils.update_frequency_list('acc_term', temp_child[0])
        #utils.update_frequency_list('don_core', temp_child[1])
        #utils.update_frequency_list('acc_core', temp_child[2])

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
        list of all possible unit SMILES with format: [acc_term_L, don_core, acc_core, acc_term_R]
    '''

    # set up data directories
    setup = subprocess.call('mkdir /ihome/ghutchison/blp62/GA/running_GA/opt_xyz /ihome/ghutchison/blp62/GA/running_GA/GFN2_input /ihome/ghutchison/blp62/GA/running_GA/GFN2_output /ihome/ghutchison/blp62/GA/running_GA/generations /ihome/ghutchison/blp62/GA/running_GA/sTDDFT_output /ihome/ghutchison/blp62/GA/running_GA/sTDDFT_input', shell=True)

    # create directory for generation 1
    gen1 = subprocess.call('mkdir /ihome/ghutchison/blp62/GA/running_GA/generations/%s' % (gen_counter), shell=True)

    # makes initial unit frequency list
    utils.init_freq_lists(smiles_list)

    # create inital population as list of NFAs
    population = []
    population_str = []
    counter = 0
    while counter < pop_size:
        # creates list to append the unit indices to
        temp_NFA = []

        # select acceptor terminal group for NFA
        acc_term_index = random.randint(0, len(smiles_list[0])-1)
        temp_NFA.append(acc_term_index)

        # select donor core group for NFA
        don_core_index = random.randint(0, len(smiles_list[1])-1)
        temp_NFA.append(don_core_index)

        # select acceptor core group for NFA
        acc_core_index = random.randint(0, len(smiles_list[2])-1)

        # making sure that the same core and terminal units won't be picked (necessary for A-D-A'-D-A structure)
        while acc_core_index == acc_term_index:
            acc_core_index = random.randint(0, len(smiles_list[2])-1)
        temp_NFA.append(acc_core_index)

        # make SMILES string of NFA
        temp_NFA_str = utils.make_NFA_str(acc_term_index, don_core_index, acc_core_index, smiles_list)

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
            # run xtb and sTD-DFT
            run_geom_opt(temp_NFA, canonical_smiles, gen_counter)

            # updates frequency lists
            utils.update_frequency_list('acc_term', temp_NFA[0])
            utils.update_frequency_list('don_core', temp_NFA[1])
            utils.update_frequency_list('acc_core', temp_NFA[2])
            counter += 1  


            '''filename = utils.make_file_name(temp_NFA)

            # checks to see if weird geometry during xtb
            if utils.check_geom_opt(canonical_smiles, '/ihome/ghutchison/blp62/GA/running_GA/opt_xyz/%s.xyz'% (filename)) == False:
                pass
            else:
                # appends the list of unit indices as [acc_term_index, don_core_index, acc_core_index]
                population.append(temp_NFA) 
                # appends the SMILES of the new NFA
                population_str.append(canonical_smiles) '''


    # create new analysis output files
    with open('/ihome/ghutchison/blp62/GA/running_GA/quick_analysis_data.csv', 'w') as quick_file:
        # write analysis file headers
        quick_file.write('gen,min_PCE, max_PCE, med_PCE, min_deltaHOMO,max_deltaHOMO,med_deltaHOMO,min_summedoscs,max_summedoscs,med_summedoscs\n')

    with open('/ihome/ghutchison/blp62/GA/running_GA/full_analysis_data.csv', 'w') as analysis_file:
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

    #path = '/ihome/ghutchison/blp62/GA/running_GA/generations/%s' %(prev_gen)
    #list_subfolders = [f.path for f in os.scandir(path) if f.is_dir()]

    #for x in list_subfolders:
        #new_path = path + '/' + x
        
        #TODO: add check for geom opt. If bad, run function to add more molecules until pop_size reached
        #utils.check_geom_opt(canonical_smiles, )
            
    #save_opt_file = subprocess.call('(cp /ihome/ghutchison/blp62/GA/running_GA/GFN2_input/%s/xtbopt.xyz /ihome/ghutchison/blp62/GA/running_GA/opt_xyz/%s.xyz)' % (file_name, file_name), shell=True)

    # copy GFN2 output file ints GFN2_output directory
    #save_GFN2_file = subprocess.call('(cp /ihome/ghutchison/blp62/GA/running_GA/GFN2_input/%s/%s.out /ihome/ghutchison/blp62/GA/running_GA/GFN2_output/%s.out)' % (file_name, file_name, file_name), shell=True)

    # delete xtb run directory for the NFA
    #del_NFAdir = subprocess.call('(rm -r /ihome/ghutchison/blp62/GA/running_GA/GFN2_input/%s)' % (file_name), shell=True)    

    # removes test.out file
    remove_test = subprocess.call('(cd /ihome/ghutchison/blp62/GA/running_GA/generations/%s && rm /ihome/ghutchison/blp62/GA/running_GA/generations/%s/test.out)' % (prev_gen, prev_gen), shell=True)
    
    for output in glob.iglob('/ihome/ghutchison/blp62/GA/running_GA/generations/%s/*.out' % (prev_gen)):
        
        NFA_filename = output.split('/')[-1].split('.')[0]

        orig_smi = utils.xyz_to_smiles('/ihome/ghutchison/blp62/GA/running_GA/GFN2_input/%s.xyz' % (NFA_filename))
        # check if geom is good. If not, do not evaluate and pass to next generation
        if utils.check_geom_opt(orig_smi, '/ihome/ghutchison/blp62/GA/running_GA/opt_xyz/%s.xyz'% (NFA_filename)) == False:
            copy_badgeom = subprocess.call('(cp /ihome/ghutchison/blp62/GA/running_GA/opt_xyz/%s.xyz /ihome/ghutchison/blp62/GA/running_GA/bad_geom_xyz/%s.xyz)' % (NFA_filename, NFA_filename), shell=True)
            pass
        else:

            prev_pop_filename.append(NFA_filename)
        
            # move a copy of the sTD-DFT output files into an overall directory
            #move_files(NFA_filename, prev_gen)

            # turn string of NFA indices into list of indices
            NFA = utils.make_NFA_from_filename(NFA_filename)
            prev_pop.append(NFA)

    # parses oscillator strength and sum of oscillator strengths from output files
    props_list = find_elec_prop(prev_pop_filename)

    deltaHOMO = props_list[0]
    sum_oscs = props_list[1]

    # find PCE and names of NFA. fitness_list = [ranked_NFA_names, ranked_PCE, ranked_deltaHOMO, ranked_summed_oscs, ranked_best_donor]
    fitness_list = scoring.PCE_prediction(deltaHOMO, sum_oscs, prev_pop)

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
 
    with open('/ihome/ghutchison/blp62/GA/running_GA/quick_analysis_data.csv', 'a') as quick_file:
        # write to quick analysis file
        quick_file.write('%d,%f,%f,%f,%f,%f,%f,%f,%f,%f\n' % (prev_gen, min_PCE, max_PCE, med_PCE, min_test_deltaHOMO, max_test_deltaHOMO, med_test_deltaHOMO, min_test_summedoscs, max_test_summedoscs, med_test_summedoscs))

    # find values and write to full analysis file
    # loop over every NFA in population
    for x in range(len(fitness_list[0])):
        NFA = fitness_list[0][x]
        file_name = utils.make_file_name(NFA)
        donor = fitness_list[4][x]
        PCE = fitness_list[1][x]
        deltaHOMO = fitness_list[2][x][0]
        summed_oscs = fitness_list[3][x]
        
        with open('/ihome/ghutchison/blp62/GA/running_GA/full_analysis_data.csv', 'a') as analysis_file:
            analysis_file.write('%d,%s,%f,%f,%f,%s,\n' % (prev_gen, file_name, deltaHOMO, summed_oscs, PCE, donor))

    # make backup copies of output files
    shutil.copy('/ihome/ghutchison/blp62/GA/running_GA/quick_analysis_data.csv', '/ihome/ghutchison/blp62/GA/running_GA/quick_analysis_data_copy.csv')
    shutil.copy('/ihome/ghutchison/blp62/GA/running_GA/full_analysis_data.csv', '/ihome/ghutchison/blp62/GA/running_GA/full_analysis_data_copy.csv')

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
        [acceptor_terminal_left, donor_core, acceptor_core, acceptor_terminal_right]
    fitness_list: list
        fitness_list = [ranked_NFA_names, ranked_PCE, ranked_deltaHOMO, ranked_summed_oscs, ranked_best_donor]
    '''

    # makes a directory in generations directory with current genetation number
    setup = subprocess.call('mkdir /ihome/ghutchison/blp62/GA/running_GA/generations/%s' % (gen_counter), shell=True)


    # load openbabel
    #load_obabel = subprocess.call('module load intel openbabel', shell=True)

    # Selection - select heaviest (best) 50% of NFAs as parents
    parents = parent_select(fitness_list, pop_size)

    new_pop = deepcopy(parents)
    new_pop_str = []
    for parent in new_pop:
        parent_str = utils.make_NFA_str(int(parent[0]), int(parent[1]), int(parent[2]), smiles_list)
        new_pop_str.append(parent_str)
        # make file name string w/ convention   
        file_name = utils.make_file_name(parent)
        duplicate_sTDDFT = subprocess.call('(cp /ihome/ghutchison/blp62/GA/running_GA/sTDDFT_output/%s.out /ihome/ghutchison/blp62/GA/running_GA/generations/%s/%s.out)' % (file_name, gen_counter, file_name), shell=True)
        duplicate_sTDDFT_prop = subprocess.call('(cp /ihome/ghutchison/blp62/GA/running_GA/sTDDFT_output/%s.prop /ihome/ghutchison/blp62/GA/running_GA/generations/%s/%s.prop)' % (file_name, gen_counter, file_name), shell=True)
        #run_sTDDFT(gen_counter, parent)


    while len(new_pop) < pop_size:

        cross_mut_pop = crossover_mutate(parents, pop_size, smiles_list, new_pop_str)
        # checks for duplication in generation
        while cross_mut_pop == False:
            cross_mut_pop = crossover_mutate(parents, pop_size, smiles_list, new_pop_str)
        
        new_NFA = cross_mut_pop[0]
        new_NFA_str = cross_mut_pop[1]

        run_geom_opt(new_NFA, new_NFA_str, gen_counter)

        new_pop.append(new_NFA)
        new_pop_str.append(new_NFA_str)

        #filename = utils.make_file_name(new_NFA)

        utils.update_frequency_list('acc_term', new_NFA[0])
        utils.update_frequency_list('don_core', new_NFA[1])
        utils.update_frequency_list('acc_core', new_NFA[2])



        # checks to see if weird geometry during xtb
        '''while utils.check_geom_opt(new_NFA_str, '/ihome/ghutchison/blp62/GA/running_GA/opt_xyz/%s.xyz'% (filename)) == False:
            cross_mut_pop = crossover_mutate(parents, pop_size, smiles_list, new_pop_str)
            while cross_mut_pop == False:
                cross_mut_pop = crossover_mutate(parents, pop_size, smiles_list, new_pop_str)
            new_NFA = cross_mut_pop[0]
            new_NFA_str = cross_mut_pop[1]
            run_geom_opt(new_NFA, new_NFA_str, gen_counter)
            filename = utils.make_file_name(new_NFA)'''
        #run_sTDDFT(gen_counter, new_NFA)

    # updates the frequency lists
    '''NFAs = []
    for output in glob.iglob('/ihome/ghutchison/blp62/GA/running_GA/generations/%s/*.inp' % gen_counter):
            
        NFA_filename = output.split('/')[-1].split('.')[0]

        if NFA_filename not in NFAs:
            NFAs.append(NFA_filename)
            units = NFA_filename.split('_')
            utils.update_frequency_list('acc_term', int(units[0]))
            utils.update_frequency_list('don_core', int(units[1]))
            utils.update_frequency_list('acc_core', int(units[2]))'''
'''
def add_new_NFA(prev_gen):
    if prev_gen == 1:
        # makes initial unit frequency list
        utils.init_freq_lists(smiles_list)

        # creates list to append the unit indices to
        temp_NFA = []

        # select acceptor terminal group for NFA
        acc_term_index = random.randint(0, len(smiles_list[0])-1)
        temp_NFA.append(acc_term_index)

        # select donor core group for NFA
        don_core_index = random.randint(0, len(smiles_list[1])-1)
        temp_NFA.append(don_core_index)

        # select acceptor core group for NFA
        acc_core_index = random.randint(0, len(smiles_list[2])-1)

        # making sure that the same core and terminal units won't be picked (necessary for A-D-A'-D-A structure)
        while acc_core_index == acc_term_index:
            acc_core_index = random.randint(0, len(smiles_list[2])-1)
        temp_NFA.append(acc_core_index)

        # make SMILES string of NFA
        temp_NFA_str = utils.make_NFA_str(acc_term_index, don_core_index, acc_core_index, smiles_list)

        try:
            # convert to canonical SMILES to check for duplication
            canonical_smiles = Chem.MolToSmiles(Chem.MolFromSmiles(temp_NFA_str))
        except:
            # prevents molecules with incorrect valence, which canonical smiles will catch and throw error
            print(temp_NFA_str)
            print('Incorrect valence, could not perform canonical smiles')
            pass

        # check for duplication
        filename = utils.make_file_name(temp_NFA)

        # TODO: fix this exists line. 
        '/ihome/ghutchison/blp62/GA/running_GA/generations/%s/' %(prev_gen) + filename + '.out'
        exists_in_gen = os.path.isfile('/ihome/ghutchison/blp62/GA/running_GA/generations/%s/' %(prev_gen) + filename + '.out')
        if filename in population_st:
            pass
        else:
            #make_xtb_files(temp_NFA, canonical_smiles, gen_counter) 
            run_geom_opt(temp_NFA, canonical_smiles, gen_counter)  

    else:


def check_calcs(path, list_subfolders, unit_list, prev_gen):
    # goes through the list of subfolders in generations
    for x in list_subfolders:
        # file path can now be '/ihome/ghutchison/blp62/GA/running_GA/generations/1'
        new_path = path + '/' + x

        # creates list of all subdirectories within generations/1
        list_running_NFAs = [f.name for f in os.scandir(new_path) if f.is_dir()]
        for i in list_running_NFAs:
            # new path as generations/1/filename/
            NFA_path = new_path + '/' + i

            # if sTD-DFT exists, add it to generation
            exists = os.path.isfile(NFA_path + '/.prop')
            if exists:

                # create SMILES from filename
                NFA = make_NFA_from_filename(i)
                NFA_SMILES = make_NFA_str(NFA[0], NFA[1], NFA[2], unit_list)

                xyz_filepath = NFA_path + '/*_xtb.xyz'
                # checks to see if geometry is good or xtb did something weird
                is_geom_good = check_geom_opt(NFA_SMILES, xyz_filepath)
                if is_geom_good == False:
                    add_new_NFA(prev_gen)
                    # TODO: remove the bad NFA so it is no longer a subdirectory
                    return False
                else:

                    # TODO: fix this area
                    duplicate_GFN2_input = subprocess.call('(cp /ihome/ghutchison/blp62/GA/running_GA/sTDDFT_output/%s.out /ihome/ghutchison/blp62/GA/running_GA/generations/%s/%s.out)' % (file_name, gen_counter, file_name), shell=True)
                    duplicate_sTDDFT_prop = subprocess.call('(cp /ihome/ghutchison/blp62/GA/running_GA/sTDDFT_output/%s.prop /ihome/ghutchison/blp62/GA/running_GA/generations/%s/%s.prop)' % (file_name, gen_counter, file_name), shell=True)
                    
                    #In this section, we want to check geom of xtb, copy all files into correct directories, rename files, and then remove directory
                    
            else:
                pass

def check_gen(gen_counter, unit_list):
    prev_gen = gen_counter - 1

    # main path to check
    path = '/ihome/ghutchison/blp62/GA/running_GA/generations/%s' % (prev_gen)
    list_subfolders = [f.name for f in os.scandir(path) if f.is_dir()]

    if len(list_subfolders) != 0:
        print('subfolders remaining. Molecule is still running cals')
        return
    else:
        check_calcs(prev_gen, path, list_subfolders)


'''

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

    # Create list of possible building block unit SMILES in specific format ([acc_term_L, don_core, acc_core, acc_term_R])
    smiles_list = utils.make_unit_list()
            
    # run if this is the initial generation
    if gen_counter == 1:
        # initialize population and submit calculations
        init_gen(pop_size, gen_counter, smiles_list)
        

        # pickle random state for restart
        randstate = random.getstate()
        rand_file = open('/ihome/ghutchison/blp62/GA/running_GA/rand_states/randstate_%s.p' % (str(gen_counter)), 'wb')
        pickle.dump(randstate, rand_file)
        rand_file.close()

    # run if NOT the initial generation
    else:
        # reload random state

        prev_gen = str(gen_counter -1)

        open_rand = open('/ihome/ghutchison/blp62/GA/running_GA/rand_states/randstate_%s.p' % (prev_gen), 'rb')
        randstate = pickle.load(open_rand)
        random.setstate(randstate)
        open_rand.close()

        # evaluate previous generation
        fitness_list = evaluate_prev_gen(gen_counter)

        # run next generation of GA
        next_gen(gen_counter, pop_size, smiles_list, fitness_list)


        # save pickle random state for restart
        randstate = random.getstate()
        rand_file = open('/ihome/ghutchison/blp62/GA/running_GA/rand_states/randstate_%s.p' % (str(gen_counter)), 'wb')
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
        for x in glob.iglob('/ihome/ghutchison/blp62/GA/running_GA/generations/%s/*.prop' % (prev_gen)):
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
    rootdir = '/ihome/ghutchison/blp62/GA/running_GA/generations/'
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
        
        
        