import numpy as np
import pandas as pd

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
