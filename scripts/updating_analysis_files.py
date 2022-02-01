import glob
import numpy as np
import pandas as pd
import csv
import shutil

# create new analysis output files
with open('/ihome/ghutchison/blp62/GA/running_GA/quick_analysis_data_UPDATED.csv', 'w') as quick_file:
    # write analysis file headers
    quick_file.write('gen,min_PCE, max_PCE, med_PCE, min_deltaHOMO,max_deltaHOMO,med_deltaHOMO,min_summedoscs,max_summedoscs,med_summedoscs\n')

with open('/ihome/ghutchison/blp62/GA/running_GA/full_analysis_data_UPDATED.csv', 'w') as analysis_file:
    # write analysis file headers
    analysis_file.write('gen,filename,deltaHOMO,summedoscs,PCE,donor\n')

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
    coreA = str(NFA[2])
    coreD = str(NFA[1])

    # make file name string for A-D-A'-D-A structure
    file_name = '%s_%s_%s_%s_%s' % (termA, coreD, coreA, coreD, termA)

    return file_name

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
        stddft_props = parse_sTDDFT(NFA)
        NFA_deltaHOMO_list.append(stddft_props[0])
        NFA_summedoscs_list.append(stddft_props[1])

    # make nested list of deltaHOMO and summed oscs lists
    stddft_prop_lists = [NFA_deltaHOMO_list, NFA_summedoscs_list]

    return stddft_prop_lists


def PCE_prediction(deltaHOMO, summed_oscs, population):
    """
    Calculates the PCE of the best acceptor-donor pair and ranks the population

    Parameters
    ----------
    deltaHOMO: list
        list of difference in energy between HOMO and HOMO-1 of each NFA in the population
    summed_oscs: list
        list of sums of oscillator strengths of each NFA in the population
    population: nested list
        nested list of unit indices that make up NFA

    Return
    ------
    ranked_population: nested list
        lists of NFAs and their PCE ranked in order of highest PCE first. Also contains deltaHOMO, summed oscs, and the best donor
        [ranked_NFA_names, ranked_PCE, ranked_deltaHOMO, ranked_summed_oscs, ranked_best_donor]
    """

    # retrieving properties from pre-calculated sTD-DFT of 20 donor polymers
    donors = pd.read_csv('/ihome/ghutchison/blp62/GA/running_GA/donor_props/donor_properties.csv')

    score_list = []
    NFA_name = []
    deltaHOMO_list = []
    summed_oscs_list = []
    best_donor = []
    
    for x in range(len(deltaHOMO)):
        # sets the acceptor properties

        deltaHOMO_acc = deltaHOMO[x][0]
        sum_oscs_acc = summed_oscs[x]

        PCE = 0
        #chooses which donor is best for the acceptor (results in highest PCE)
        for index, row in donors.iterrows():
            temp_donor_name = row['Donor']
            # equation to predict PCE
            temp_PCE = -33.08 + (1.377*sum_oscs_acc) +(4.255*deltaHOMO_acc) + (-0.4587*row['summedoscs']) + (0.1735*row['absFOM']) + (2.449*row['Electrodonating']) + (0.0009508*row['optbg'])
            
            if temp_PCE > PCE:
                PCE = temp_PCE
                donor = temp_donor_name

        score_list.append(PCE)
        NFA_name.append(population[x])
        deltaHOMO_list.append(deltaHOMO_acc)
        summed_oscs_list.append(sum_oscs_acc)
        best_donor.append(donor)

    ranked_PCE = []
    ranked_NFA_names = []
    ranked_deltaHOMO = []
    ranked_summed_oscs = []
    ranked_best_donor = []

    # make list of indicies of NFAs in population, sorted based on PCE
    ranked_indices = list(np.argsort(score_list))
    # reverse list so highest property value = 0th
    ranked_indices.reverse()

    # makes list of each property in order of highest PCE pair to lowest
    for x in ranked_indices:
        ranked_PCE.append(score_list[x])
        ranked_NFA_names.append(NFA_name[x])
        ranked_deltaHOMO.append(deltaHOMO[x])
        ranked_summed_oscs.append(summed_oscs[x])
        ranked_best_donor.append(best_donor[x])

    ranked_population = [ranked_NFA_names, ranked_PCE, ranked_deltaHOMO, ranked_summed_oscs, ranked_best_donor]

    return ranked_population


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
    # creates ordered list of unit index [acc_term, acc_core, don_core]
    NFA = [NFA_indices[0], NFA_indices[1], NFA_indices[2]]

    return NFA

def evaluate_gen(gen_counter):
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

    prev_pop = []
    prev_pop_filename = []

    print(gen_counter)

    for output in glob.iglob('/ihome/ghutchison/blp62/GA/running_GA/generations/%s/*.out' % (gen_counter)):
        print(output)
        NFA_filename = output.split('/')[-1].split('.')[0]
        prev_pop_filename.append(NFA_filename)

        print(NFA_filename)
        # turn string of NFA indices into list of indices
        NFA = make_NFA_from_filename(NFA_filename)
        prev_pop.append(NFA)

    # parses oscillator strength and sum of oscillator strengths from output files
    props_list = find_elec_prop(prev_pop_filename)

    deltaHOMO = props_list[0]
    sum_oscs = props_list[1]

    # find PCE and names of NFA. fitness_list = [ranked_NFA_names, ranked_PCE, ranked_deltaHOMO, ranked_summed_oscs, ranked_best_donor]
    fitness_list = PCE_prediction(deltaHOMO, sum_oscs, prev_pop)
    print(fitness_list)
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
 
    with open('/ihome/ghutchison/blp62/GA/running_GA/quick_analysis_data_UPDATED.csv', 'a') as quick_file:
        # write to quick analysis file
        quick_file.write('%d,%f,%f,%f,%f,%f,%f,%f,%f,%f\n' % (gen_counter, min_PCE, max_PCE, med_PCE, min_test_deltaHOMO, max_test_deltaHOMO, med_test_deltaHOMO, min_test_summedoscs, max_test_summedoscs, med_test_summedoscs))

    # find values and write to full analysis file
    # loop over every NFA in population
    for x in range(len(fitness_list[0])):
        NFA = fitness_list[0][x]
        file_name = make_file_name(NFA)
        donor = fitness_list[4][x]
        PCE = fitness_list[1][x]
        deltaHOMO = fitness_list[2][x][0]
        summed_oscs = fitness_list[3][x]
        
        with open('/ihome/ghutchison/blp62/GA/running_GA/full_analysis_data.csv_UPDATED', 'a') as analysis_file:
            analysis_file.write('%d,%s,%f,%f,%f,%s,\n' % (gen_counter, file_name, deltaHOMO, summed_oscs, PCE, donor))

    # make backup copies of output files
    shutil.copy('/ihome/ghutchison/blp62/GA/running_GA/quick_analysis_data_UPDATED.csv', '/ihome/ghutchison/blp62/GA/running_GA/quick_analysis_data_copy_UPDATED.csv')
    shutil.copy('/ihome/ghutchison/blp62/GA/running_GA/full_analysis_data_UPDATED.csv', '/ihome/ghutchison/blp62/GA/running_GA/full_analysis_data_copy_UPDATED.csv')

    # fitness_list = [ranked_NFA_names, ranked_PCE, ranked_deltaHOMO, ranked_summed_oscs, ranked_best_donor]
    return fitness_list

for gen in range(1, 15):
    evaluate_gen(gen)