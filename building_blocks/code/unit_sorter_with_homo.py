# To add a new cell, type '# %%'
# To add a new markdown cell, type '# %% [markdown]'
# %%
import cclib
import glob
import matplotlib.pyplot as plt
import csv
import os
import pandas as pd
from matplotlib import cm


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

for filename in glob.iglob('../optimized_B3LYP/*.out'):

    myfile = cclib.io.ccopen(filename)
    values = myfile.parse()

    homo = values.homos[0] # index of HOMO. assuming spin-restricted.
    homo_energy = values.moenergies[0][homo] 

    lumo = homo + 1 # index of LUMO
    lumo_energy = values.moenergies[0][lumo]

    filename = filename.split('.',1)[0]
    filename = filename.split('/')[-1]

    for smilename in glob.iglob('../smiles/*.smi'):
        smile = smilename.split('.',1)[0]
        smile = smile.split('/')[-1]
            
        if smile == filename:
            with open(smilename, "r") as smiles:
                smiles = smiles.read()
            break
    
    eng_diff = thiophene_HOMO - homo_energy
    energy_difference.append(eng_diff)

    if eng_diff < 0:
        if "t" in filename:
            donor_left_terminal_old_name.append(filename)
            index = donor_left_terminal_old_name.index(filename)
            new_name = str(index) + "dtl"
            donor_left_terminal_new_name.append(new_name)
            donor_left_terminal_smi.append(smiles)
            donor_left_terminal_homo.append(homo_energy)
            donor_left_terminal_lumo.append(lumo_energy)

            identifier = filename[0:-1]

            for right in glob.iglob('../terminal_right_smiles/*.smi'):
                right_str = right.split('.',1)[0]
                right_str = right_str.split('/')[-1]

                if identifier in right:
                    with open(right, "r") as right_smiles:
                        right_smiles = right_smiles.read()
                    donor_right_terminal_smi.append(right_smiles)
                    donor_right_terminal_old_name.append(right_str)
                    index = donor_right_terminal_old_name.index(right_str)
                    new_name = str(index) + "dtr"
                    donor_right_terminal_new_name.append(new_name)
                    donor_right_terminal_homo.append(homo_energy)
                    donor_right_terminal_lumo.append(lumo_energy)
                    break

        else:
            donor_core_smi.append(smiles)
            donor_core_old_name.append(filename)
            index = donor_core_old_name.index(filename)
            new_name = str(index) + "dc"
            donor_core_new_name.append(new_name)
            donor_core_homo.append(homo_energy)
            donor_core_lumo.append(lumo_energy)

            donor_left_terminal_smi.append(smiles)
            donor_left_terminal_old_name.append(filename)
            donor_left_terminal_new_name.append(new_name)
            donor_left_terminal_homo.append(homo_energy)
            donor_left_terminal_lumo.append(lumo_energy)

            donor_right_terminal_smi.append(smiles)
            donor_right_terminal_old_name.append(filename)
            donor_right_terminal_new_name.append(new_name)
            donor_right_terminal_homo.append(homo_energy)
            donor_right_terminal_lumo.append(lumo_energy)


    else:
        if "t" in filename:
            acceptor_left_terminal_old_name.append(filename)
            index = acceptor_left_terminal_old_name.index(filename)
            new_name = str(index) + "atl"
            acceptor_left_terminal_new_name.append(new_name)
            acceptor_left_terminal_smi.append(smiles)
            acceptor_left_terminal_homo.append(homo_energy)
            acceptor_left_terminal_lumo.append(lumo_energy)

            identifier = filename[0:-1]

            for right in glob.iglob('../terminal_right_smiles/*.smi'):
                right_str = right.split('.',1)[0]
                right_str = right_str.split('/')[-1]

                if identifier in right:
                    with open(right, "r") as right_smiles:
                        right_smiles = right_smiles.read()
                    acceptor_right_terminal_smi.append(right_smiles)
                    acceptor_right_terminal_old_name.append(right_str)
                    index = acceptor_right_terminal_old_name.index(right_str)
                    new_name = str(index) + "atr"
                    acceptor_right_terminal_new_name.append(new_name)
                    acceptor_right_terminal_homo.append(homo_energy)
                    acceptor_right_terminal_lumo.append(lumo_energy)
                    break

        else:
            acceptor_core_smi.append(smiles)
            acceptor_core_old_name.append(filename)
            index = acceptor_core_old_name.index(filename)
            new_name = str(index) + "ac"
            acceptor_core_new_name.append(new_name)
            acceptor_core_homo.append(homo_energy)
            acceptor_core_lumo.append(lumo_energy)

            acceptor_left_terminal_smi.append(smiles)
            acceptor_left_terminal_old_name.append(filename)
            acceptor_left_terminal_new_name.append(new_name)
            acceptor_left_terminal_homo.append(homo_energy)
            acceptor_left_terminal_lumo.append(lumo_energy)

            acceptor_right_terminal_smi.append(smiles)
            acceptor_right_terminal_old_name.append(filename)
            acceptor_right_terminal_new_name.append(new_name)
            acceptor_right_terminal_homo.append(homo_energy)
            acceptor_right_terminal_lumo.append(lumo_energy)


df_acc_core = pd.DataFrame({'SMILES': acceptor_core_smi, 'new name':acceptor_core_new_name, 'old name': acceptor_core_old_name, 'homo': acceptor_core_homo, 'lumo': acceptor_core_lumo})
df_acc_core.to_csv('acc_core_with_homo.csv')

df_acc_left_term = pd.DataFrame({'SMILES': acceptor_left_terminal_smi, 'new name':acceptor_left_terminal_new_name, 'old name': acceptor_left_terminal_old_name, 'homo': acceptor_left_terminal_homo, 'lumo': acceptor_left_terminal_lumo})
df_acc_left_term.to_csv('acc_left_term_with_homo.csv')

df_acc_right_term = pd.DataFrame({'SMILES': acceptor_right_terminal_smi, 'new name':acceptor_right_terminal_new_name, 'old name': acceptor_right_terminal_old_name, 'homo': acceptor_right_terminal_homo, 'lumo': acceptor_right_terminal_lumo})
df_acc_right_term.to_csv('acc_right_term_with_homo.csv')

df_don_core = pd.DataFrame({'SMILES': donor_core_smi, 'new name':donor_core_new_name, 'old name': donor_core_old_name, 'homo': donor_core_homo, 'lumo': donor_core_lumo})
df_don_core.to_csv('don_core_with_homo.csv')

df_don_left_term = pd.DataFrame({'SMILES': donor_left_terminal_smi, 'new name':donor_left_terminal_new_name, 'old name': donor_left_terminal_old_name, 'homo': donor_left_terminal_homo, 'lumo': donor_left_terminal_lumo})
df_don_left_term.to_csv('don_left_term_with_homo.csv')

df_don_right_term = pd.DataFrame({'SMILES': donor_right_terminal_smi, 'new name':donor_right_terminal_new_name, 'old name': donor_right_terminal_old_name, 'homo': donor_right_terminal_homo, 'lumo': donor_right_terminal_lumo})
df_don_right_term.to_csv('don_right_term_with_homo.csv')
   
'''
If HOMO is less than thiophene HOMO, it is donor. 
If HOMO is greater than thiphene HOMO, it is acceptor
'''



