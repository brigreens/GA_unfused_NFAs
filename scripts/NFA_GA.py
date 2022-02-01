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
#import pybel
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

def run_geom_opt(NFA, NFA_str, generation_counter):
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

    # if optimized xyz exists, skip xTB
    exists = os.path.isfile('opt_xyz/%s.xyz' % (file_name))
    if exists:
        print("xyz file existed")
        return

    # make NFA string into pybel molecule object
    mol = pybel.readstring('smi', NFA_str)
    utils.make3D(mol)


    # write NFA .xyz file to containing folder
    mol.write('xyz', 'GFN2_input/%s.xyz' % (file_name), overwrite=True)

    # make directory to run xtb in for the NFA
    mkdir_NFA = subprocess.call('(mkdir GFN2_input/%s)' % (file_name), shell=True)

    # run xTB geometry optimization
    xtb = subprocess.call('(cd GFN2_input/%s && /ihome/ghutchison/geoffh/xtb/xtb ../%s.xyz --opt > %s.out)' %
                          (file_name, file_name, file_name), shell=True)

    save_opt_file = subprocess.call(
        '(cp GFN2_input/%s/xtbopt.xyz opt_xyz/%s.xyz)' % (file_name, file_name), shell=True)

    # copy GFN2 output file ints GFN2_output directory
    save_GFN2_file = subprocess.call(
        '(cp GFN2_input/%s/%s.out GFN2_output/%s.out)' % (file_name, file_name, file_name), shell=True)

    # delete xtb run directory for the NFA
    del_NFAdir = subprocess.call('(rm -r GFN2_input/%s)' % (file_name), shell=True)    

def run_sTDDFT(gen_counter, NFA):
    '''
    Submits NFAs for sTD-DFT calculation with slurm
    Parameters
    ---------
    gen_counter: int
        current generation
    NFA: list (specific format)
        [acc_term_index, don_core_index, acc_core_index]
    '''

    # make file name string w/ convention   
    file_name = utils.make_file_name(NFA)

    # if sTD-DFT exists, add it to generation
    exists = os.path.isfile('sTDDFT_output/%s.out' % (file_name))
    if exists:
        duplicate_sTDDFT = subprocess.call('(cp sTDDFT_output/%s.out generations/%s/%s.out)' % (file_name, gen_counter, file_name), shell=True)
        duplicate_sTDDFT_prop = subprocess.call('(cp sTDDFT_output/%s.prop generations/%s/%s.prop)' % (file_name, gen_counter, file_name), shell=True)

        print("sTDDFT file existed. Copied to current generation")
        return 
    else:
        # move optimized xyz file into sTDDFT input folder
        create_input = subprocess.call('(cd generations/%s && obabel ../../opt_xyz/%s.xyz -O .inp -o orcainp -m -xf ../../cam_b3lyp_def2_SVP_pp.txt)' % (gen_counter, file_name), shell=True)
        # submit sTD-DFT calculation
        sTDDFT_submit = subprocess.call('(cd generations/%s && sbatch -J %s.inp ~blp62/orca4_sTDDFT_pp.slurm)' % (gen_counter, file_name), shell=True)

def move_files(file_name, prev_gen):
    '''
    moves output files from previous generation into organized directories

    Parameters
    ----------
    file_name: str
        NFA file name (w/o extension) showing unit indicies
        e.g. 100_200_230_200_100 for a certain NFA with A-D-A'-D-A structure
    '''

    # makes a copy of the sTDDFT output and input files into separate directories
    copy_sTDDFT_output = subprocess.call('(cp generations/%s/%s.out sTDDFT_output/%s.out)' % (prev_gen, file_name, file_name), shell=True)
    copy_sTDDFT_prop = subprocess.call('(cp generations/%s/%s.prop sTDDFT_output/%s.prop)' % (prev_gen, file_name, file_name), shell=True)

    #copy_sTDDFT_input = subprocess.call('(cp %s.inp sTDDFT_input/%s.inp)' % (file_name, file_name), shell=True)


def find_elec_prop(file_name):
    '''
    Calculates oscillator strength and sum of oscillator strength of each NFA in population

    Parameters
    ---------
    population: list

    Returns
    -------
    calc_prop_lists: list
        nested list of [list of oscillator strength, list of sum of oscillator strengths]
    '''

    NFA_deltaHOMO_list = []
    NFA_summedoscs_list = []

    # parse sTD-DFT output files
    for NFA in file_name:
        # parses the sTD-DFT output file to find the first oscillator strength and sum of all oscillator strengths 
        stddft_props = utils.parse_sTDDFT(NFA)
        NFA_deltaHOMO_list.append(stddft_props[0])
        NFA_summedoscs_list.append(stddft_props[1])

    # make nested list of deltaHOMO and summed oscs lists
    stddft_prop_lists = [NFA_deltaHOMO_list, NFA_summedoscs_list]

    return stddft_prop_lists


def parent_select(fitness_list):
    '''
    Finds top half of population

    Parameters
    ---------
    population: list
        list of NFAs in population
    fitness_list: list
        fitness_list = [ranked_NFA_names, ranked_PCE, ranked_deltaHOMO, ranked_summed_oscs, ranked_best_donor]

    Returns
    -------
    parent_list: list
        top half of population
    '''
    # find number of parents (half of population)
    parent_count = int(len(fitness_list[0]) / 2)

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
    else:
        # updates the unit frequency list for the new NFAs
        utils.update_frequency_list('acc_term', temp_child[0])
        utils.update_frequency_list('don_core', temp_child[1])
        utils.update_frequency_list('acc_core', temp_child[2])

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
    setup = subprocess.call('mkdir opt_xyz GFN2_input GFN2_output generations sTDDFT_output sTDDFT_input', shell=True)

    gen1 = subprocess.call('mkdir generations/%s' % (gen_counter), shell=True)
            
    # load openbabel
    load_obabel = subprocess.call('module load intel openbabel', shell=True)


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
        # increase frequency count for unit
        utils.update_frequency_list('acc_term', acc_term_index )
        
        # select donor core group for NFA
        don_core_index = random.randint(0, len(smiles_list[1])-1)
        # increase frequency count for unit in don_core_list
        utils.update_frequency_list('don_core', don_core_index )
        temp_NFA.append(don_core_index)

        # select acceptor core group for NFA
        acc_core_index = random.randint(0, len(smiles_list[2])-1)

        # making sure that the same core and terminal units won't be picked (necessary for A-D-A'-D-A structure)
        while acc_core_index == acc_term_index:
            acc_core_index = random.randint(0, len(smiles_list[2])-1)
        # increase frequency count for unit in acc_core_list
        utils.update_frequency_list('acc_core', acc_core_index )
        temp_NFA.append(acc_core_index)

        # make SMILES string of NFA
        temp_NFA_str = utils.make_NFA_str(acc_term_index, don_core_index, acc_core_index, smiles_list)

        print(temp_NFA)

        try:
            # convert to canonical SMILES to check for duplication
            canonical_smiles = Chem.MolToSmiles(Chem.MolFromSmiles(temp_NFA_str))
        except:
            pass

        # check for duplication
        if canonical_smiles in population_str:
            pass
        else:
            run_geom_opt(temp_NFA, canonical_smiles, gen_counter)  

            filename = utils.make_file_name(temp_NFA)

            # checks to see if weird geometry during xtb
            if utils.check_geom_opt(canonical_smiles, 'opt_xyz/%s.xyz'% (filename)) == False:
                pass
            else:
                # appends the list of unit indices as [acc_term_index, don_core_index, acc_core_index]
                population.append(temp_NFA) 
                # appends the SMILES of the new NFA
                population_str.append(canonical_smiles) 

                run_sTDDFT(gen_counter, temp_NFA)

                counter += 1


    # create new analysis output files
    with open('quick_analysis_data.csv', 'w') as quick_file:
        # write analysis file headers
        quick_file.write('gen,min_PCE, max_PCE, med_PCE, min_deltaHOMO,max_deltaHOMO,med_deltaHOMO,min_summedoscs,max_summedoscs,med_summedoscs\n')

    with open('full_analysis_data.csv', 'w') as analysis_file:
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
    remove_test = subprocess.call('(cd generations/%s && rm test.out)' % (prev_gen), shell=True)
    for output in glob.iglob('generations/%s/*.out' % (prev_gen)):
        
        NFA_filename = output.split('/')[-1].split('.')[0]
        prev_pop_filename.append(NFA_filename)
    
        # move a copy of the sTD-DFT output files into an overall directory
        move_files(NFA_filename, prev_gen)

        # turn string of NFA indices into list of indices
        NFA = utils.make_NFA_from_filename(NFA_filename)
        prev_pop.append(NFA)

    # parses oscillator strength and sum of oscillator strengths from output files
    props_list = find_elec_prop(prev_pop_filename)

    deltaHOMO = props_list[0]
    sum_oscs = props_list[1]

    #switch_dir = subprocess.call('cd ..', shell=True)

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
 
    with open('quick_analysis_data.csv', 'a') as quick_file:
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
        
        with open('full_analysis_data.csv', 'a') as analysis_file:
            analysis_file.write('%d,%s,%f,%f,%f,%s,\n' % (prev_gen, file_name, deltaHOMO, summed_oscs, PCE, donor))

    # make backup copies of output files
    shutil.copy('quick_analysis_data.csv', 'quick_analysis_data_copy.csv')
    shutil.copy('full_analysis_data.csv', 'full_analysis_data_copy.csv')

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
    setup = subprocess.call('mkdir generations/%s' % (gen_counter), shell=True)

    # load openbabel
    load_obabel = subprocess.call('module load intel openbabel', shell=True)

    # Selection - select heaviest (best) 50% of NFAs as parents
    parents = parent_select(fitness_list)

    new_pop = deepcopy(parents)
    new_pop_str = []
    for parent in new_pop:
        parent_str = utils.make_NFA_str(int(parent[0]), int(parent[1]), int(parent[2]), smiles_list)
        new_pop_str.append(parent_str)
        run_sTDDFT(gen_counter, parent)


    while len(new_pop) < pop_size:

        cross_mut_pop = crossover_mutate(parents, pop_size, smiles_list, new_pop_str)
        # checks for duplication in generation
        while cross_mut_pop == False:
            cross_mut_pop = crossover_mutate(parents, pop_size, smiles_list, new_pop_str)
        
        new_NFA = cross_mut_pop[0]
        new_NFA_str = cross_mut_pop[1]

        run_geom_opt(new_NFA, new_NFA_str, gen_counter)

        filename = utils.make_file_name(new_NFA)

        # checks to see if weird geometry during xtb
        while utils.check_geom_opt(new_NFA_str, 'opt_xyz/%s.xyz'% (filename)) == False:
            cross_mut_pop = crossover_mutate(parents, pop_size, smiles_list, new_pop_str)
            while cross_mut_pop == False:
                cross_mut_pop = crossover_mutate(parents, pop_size, smiles_list, new_pop_str)
            new_NFA = cross_mut_pop[0]
            new_NFA_str = cross_mut_pop[1]
            run_geom_opt(new_NFA, new_NFA_str, gen_counter)
            filename = utils.make_file_name(new_NFA)

        new_pop.append(new_NFA)
        new_pop_str.append(new_NFA_str)
        run_sTDDFT(gen_counter, new_NFA)


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
        rand_file = open('rand_states/randstate_%s.p' % (str(gen_counter)), 'wb')
        pickle.dump(randstate, rand_file)
        rand_file.close()

    # run if NOT the initial generation
    else:
        # reload random state

        prev_gen = str(gen_counter -1)

        open_rand = open('rand_states/randstate_%s.p' % (prev_gen), 'rb')
        randstate = pickle.load(open_rand)
        random.setstate(randstate)
        open_rand.close()

        # evaluate previous generation
        fitness_list = evaluate_prev_gen(gen_counter)

        # run next generation of GA
        next_gen(gen_counter, pop_size, smiles_list, fitness_list)


        # save pickle random state for restart
        randstate = random.getstate()
        rand_file = open('rand_states/randstate_%s.p' % (str(gen_counter)), 'wb')
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

    if gen_counter == 1:
        return True
    else:
        prev_gen = str(gen_counter -1)
        pop_size = 32

        # checks to see if all of the sTDDFT calculations are finished. *.prop used since it only appears when calculation is finished. 
        prop_counter = 0
        for x in glob.iglob('generations/%s/*.prop' % (prev_gen)):
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

    # gen_counter is the generation you would like to run next 
    gen_counter = args.gen

    # checks to make sure all of the previous generation's sTD-DFT files are complete
    if check_prev_gen_complete(gen_counter) == True:
        main(gen_counter)
    else:
        print('Previous generation has not finished sTD-DFT calculations yet')