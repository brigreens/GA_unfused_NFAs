#!/usr/bin/env python
# coding: utf-8

from __future__ import print_function
import pandas as pd
import numpy as np
import glob
from openbabel import pybel
import scipy.constants as constants
from rdkit import Chem
import os

def parse_sTDDFT(filename):
    '''
    Parses through sTD-DFT output files

    Parameters
    -----------
    filename: str
        path to output file
    donor_or_acc: str
        define if molecule is acceptor ('acc') or donor ('don')

    Returns
    -------
    acceptors_sTD[filename] or donors_sTD[filename]: list
        adds list of sTDDFT descriptors to a dictionary for the molecule
    '''
    outputs = []
    with open(filename, 'r', encoding = 'utf-8') as file:
        line = file.readline()
        oscs = []
        wavelength = []
        energyEV = []
        wavenumber = []
        while line:
            if 'ordered frontier orbitals' in line:
                for x in range(11):
                    line = file.readline()
                HOMOminus1 = float(line[9:15])
                line = file.readline()
                HOMO = float(line[9:15])
                line = file.readline()
                line = file.readline()
                LUMO = float(line[9:15])
                line = file.readline()
                LUMOplus1 = float(line[9:15])

                deltaHOMO = abs(HOMOminus1 - HOMO)
                deltaLUMO = abs(LUMO - LUMOplus1)

            elif 'ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS' in line:
                for x in range(5):
                    line = file.readline()
                opt_bg = float(line[6:15])                    
                while line != '\n':
                    oscs.append(float(line[25:37]))
                    wavelength.append(float(line[17:23]))
                    energyEV.append(float(line[6:15]) * 1.2398e-4)
                    wavenumber.append(float(line[6:15]))
                    line = file.readline()
                    
            elif 'FINAL SINGLE POINT ENERGY' in line:
                SinglePointEn = float(line[-22:-1])

            elif 'Magnitude (Debye)' in line:                    
                dipmom = float(line[-9:-1])

            line = file.readline()  
        line = file.readline()
   
    if len(oscs) != 0:

        #Finds the lowest energy strong oscillator strength (Osc >= 0.1) and its corresponding energy    
        for i in range(len(oscs)+1): 
            if i < len(oscs): 
                if oscs[i] >= 0.1:
                    index_of_oscs = i
                    break
                else:
                    continue
            else: #if the spectrum has no oscillation strengths greater than 0.1
                maxoscs = max(oscs) #finds the maximum osc. strength
                index_of_oscs = oscs.index(maxoscs)
        transition_energy = wavenumber[index_of_oscs] #transition energy in units of cm^-1 of of lowest energy transiton with largest osc. strength
        strongest_osc = oscs[index_of_oscs]  #lowest energy osillator strength greater than 0.1 (or max osc. strength is all less than 0.1)

        highest_oscs = 0.0
        first_oscs = oscs[0]
        firstenergytransitioneV = energyEV[0]
        firstenergytransitionwavenumber = wavenumber[0]
        if len(oscs) < 3:
            for i in range(len(oscs)):
                if  oscs[i] > highest_oscs:
                    highest_oscs = oscs[i]
                    lowestenergytransitioneV= energyEV[i]
                    lowestenergytransitionwavenumber= wavenumber[i]

        else:
            for x in range(3):
                if  oscs[x] > highest_oscs:
                    highest_oscs = oscs[x]
                    lowestenergytransitioneV= energyEV[x]
                    lowestenergytransitionwavenumber= wavenumber[x]


        summed_oscs = np.sum(oscs)

        PCE = -33.08 + (1.377*summed_oscs ) +(4.255*deltaHOMO) + (-0.4587*14.24515589) + (0.1735*62.75526668) + (2.449*6.2001) + (0.0009508*17569.8)

        outputs.extend((HOMOminus1, HOMO, LUMO, LUMOplus1, deltaHOMO, deltaLUMO, opt_bg, strongest_osc, SinglePointEn, dipmom, summed_oscs, first_oscs, highest_oscs, firstenergytransitioneV, firstenergytransitionwavenumber, lowestenergytransitioneV, lowestenergytransitionwavenumber, PCE))

        filename = filename.split('/')[-1].split('.',1)[0]

        acceptors_sTD[filename] = outputs  

        return acceptors_sTD[filename]


acceptors_sTD = {}

all_files = []
rootdir = '/ihome/ghutchison/blp62/GA/running_GA/opt_GA/generations'
for subdir, dirs, files in os.walk(rootdir):
    for file in files:
        if file.endswith(".out"):
            if file not in all_files:
                filename = os.path.join(subdir, file)
                print(filename)
                parse_sTDDFT(filename)
                all_files.append(file)

    
df_acc = pd.DataFrame.from_dict(acceptors_sTD, orient = 'index', columns = ['HOMO-1 (eV)', 'HOMO (eV)', 'LUMO (eV)', 'LUMO+1 (eV)', 'Delta HOMO', 'delta LUMO','optical bandgap (cm-1)', 'oscillator strength', 'single point energy', 'dipole moment (debye)', 'summed oscs', 'first oscs', 'highest oscs under ten', 'first Energy Transition eV', 'first Energy transition wavenumber', 'lowest Energy Transition eV', 'lowest Energy transition wavenumber', 'PCE'])
df_acc = df_acc.rename_axis('Molecule')
df_acc.to_csv('/ihome/ghutchison/blp62/GA/running_GA/opt_GA/GA_data_v4.csv')  




