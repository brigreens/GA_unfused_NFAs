import cclib
import glob
import matplotlib.pyplot as plt
import csv
import os
import pandas as pd


def append_list_as_row(file_name, list_of_elem):
    # Open file in append mode
    #with open(file_name, 'a+', newline='') as write_obj:
    with open(file_name, 'a+', newline='') as write_obj: 
        # Create a writer object from csv module
        csv_writer = csv.writer(write_obj, quoting = csv.QUOTE_NONE, escapechar=',')
        # Add contents of list as last row in the csv file
        csv_writer.writerow(list_of_elem)

thiophene_HOMO = -6.5467

acceptor_core_smi = []
acceptor_core_old_name = []
acceptor_core_new_name = []
acceptor_core_homo = []
acceptor_core_lumo = []

acceptor_left_terminal_smi = []
acceptor_left_terminal_old_name = []
acceptor_left_terminal_new_name = []
acceptor_left_terminal_homo = []
acceptor_left_terminal_lumo = []

acceptor_right_terminal_smi = []
acceptor_right_terminal_old_name = []
acceptor_right_terminal_new_name = []
acceptor_right_terminal_homo = []
acceptor_right_terminal_lumo = []

donor_core_smi = []
donor_core_old_name = []
donor_core_new_name = []
donor_core_homo = []
donor_core_lumo = []

donor_left_terminal_smi = []
donor_left_terminal_old_name = []
donor_left_terminal_new_name = []
donor_left_terminal_homo = []
donor_left_terminal_lumo = []

donor_right_terminal_smi = []
donor_right_terminal_old_name = []
donor_right_terminal_new_name = []
donor_right_terminal_homo = []
donor_right_terminal_lumo = []

energy_difference = []

don_core_count = 62 # index of last donor core unit
don_term_count = 87 # index of last donor terminal unit
acc_core_count = 42 # index of last acceptor core unit
acc_term_count = 102 # index of last acceptor terminal unit

for filename in glob.iglob('/ihome/ghutchison/blp62/GA/running_GA/more_units/new_units/*.out'):

    myfile = cclib.io.ccopen(filename)
    values = myfile.parse()

    homo = values.homos[0] # index of HOMO. assuming spin-restricted.
    homo_energy = values.moenergies[0][homo] 

    lumo = homo + 1 # index of LUMO
    lumo_energy = values.moenergies[0][lumo]

    bandgap = lumo_energy-homo_energy

    eng_diff = thiophene_HOMO - homo_energy
    energy_difference.append(eng_diff)

    filename = filename.split('.',1)[0]
    filename = filename.split('/')[-1]

    with open('/ihome/ghutchison/blp62/GA/running_GA/more_units/unit_csvs/1235MonomerList.csv', newline='') as new_mons:
        count = 0
        for smiles in new_mons:
            if str(count) == filename:
                unit_smiles = smiles
                break
            count += 1


    if eng_diff < 0:
        don_core_count +=1
        don_term_count +=1
        append_list_as_row('/ihome/ghutchison/blp62/GA/running_GA/more_units/unit_csvs/don_core_with_homo.csv', [don_core_count, unit_smiles, don_core_count, filename, homo_energy, lumo_energy, bandgap])
        append_list_as_row('/ihome/ghutchison/blp62/GA/running_GA/more_units/unit_csvs/don_left_term_with_homo.csv', [don_term_count, unit_smiles, don_term_count, filename, homo_energy, lumo_energy, bandgap])
        append_list_as_row('/ihome/ghutchison/blp62/GA/running_GA/more_units/unit_csvs/don_right_term_with_homo.csv', [don_term_count, unit_smiles, don_term_count, filename, homo_energy, lumo_energy, bandgap])
        

    else:
        acc_core_count +=1
        acc_term_count +=1
        append_list_as_row('/ihome/ghutchison/blp62/GA/running_GA/more_units/unit_csvs/acc_core_with_homo.csv', [acc_core_count, unit_smiles, acc_core_count, filename, homo_energy, lumo_energy, bandgap])
        append_list_as_row('/ihome/ghutchison/blp62/GA/running_GA/more_units/unit_csvs/acc_left_term_with_homo.csv', [acc_term_count, unit_smiles, acc_term_count, filename, homo_energy, lumo_energy, bandgap])
        append_list_as_row('/ihome/ghutchison/blp62/GA/running_GA/more_units/unit_csvs/acc_right_term_with_homo.csv', [acc_term_count, unit_smiles, acc_term_count, filename, homo_energy, lumo_energy, bandgap])