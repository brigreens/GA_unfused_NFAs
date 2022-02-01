import numpy as np
import pandas as pd
import glob
import csv

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

    with open(file_name, 'r', encoding = 'utf-8') as file:
        line = file.readline()
        oscs = []
        
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

                deltaHOMO = (abs(HOMOminus1 - HOMO))

            if 'ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS' in line:
                for x in range(5):
                    line = file.readline()                
                while line != '\n':
                    oscs.append(float(line[25:37]))
                    line = file.readline()
            line = file.readline()  
        line = file.readline()
   
    summed_oscs = np.sum(oscs)

    stddft_prop = [deltaHOMO, summed_oscs, HOMO, LUMO]

    return stddft_prop

def append_list_as_row(file_name, list_of_elem):
    # Open file in append mode
    with open(file_name, 'a+', newline='') as write_obj: 
        # Create a writer object from csv module
        csv_writer = csv.writer(write_obj, quoting = csv.QUOTE_NONE, escapechar=',')
        # Add contents of list as last row in the csv file
        csv_writer.writerow(list_of_elem)

# retrieving properties from pre-calculated sTD-DFT of 20 donor polymers
donors = pd.read_csv('/ihome/ghutchison/blp62/GA/running_GA/donor_props/donor_properties.csv')
path = '/ihome/ghutchison/blp62/GA/running_GA/ADDDA/PCE_predictions.csv'


for output in glob.iglob('/ihome/ghutchison/blp62/GA/running_GA/ADDDA/output_files/*.out'):
        
    outputs = parse_sTDDFT(output)
    delta_HOMO = outputs[0]
    summed_oscs = outputs[1]
    HOMO = outputs[2]
    LUMO = outputs[3]

    NFA_filename = output.split('/')[-1].split('.')[0]

    for index, row in donors.iterrows():
        temp_donor_name = row['Donor']
        # equation to predict PCE
        temp_PCE = -33.08 + (1.377*summed_oscs ) +(4.255*delta_HOMO) + (-0.4587*row['summedoscs']) + (0.1735*row['absFOM']) + (2.449*row['Electrodonating']) + (0.0009508*row['optbg'])
        
        append_list_as_row(path, [NFA_filename, temp_donor_name, delta_HOMO, summed_oscs, HOMO, LUMO, temp_PCE])
