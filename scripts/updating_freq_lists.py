import glob
import csv
import pandas as pd


def make_unit_list():
    '''
    Makes a list of lists of the SMILES of different types of building block units

    Returns
    -------
    unit_list: list
        list of lists of SMILES
        [acceptor_terminal_left, donor_core, acceptor_core, acceptor_terminal_right]
    '''

    acc_term_L = []
    don_core = []
    acc_core = []
    acc_term_R = []

    with open("/ihome/ghutchison/blp62/GA/running_GA/building_blocks/acc_core.csv", "r") as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        next(csv_reader)
        for row in csv_reader:
            acc_core.append(row[1])

    with open("/ihome/ghutchison/blp62/GA/running_GA/building_blocks/acc_left_term.csv", "r") as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        next(csv_reader)
        for row in csv_reader:
            acc_term_L.append(row[1])

    with open("/ihome/ghutchison/blp62/GA/running_GA/building_blocks/acc_right_term.csv", "r") as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        next(csv_reader)
        for row in csv_reader:
            acc_term_R.append(row[1])

    with open("/ihome/ghutchison/blp62/GA/running_GA/building_blocks/don_core.csv", "r") as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        next(csv_reader)
        for row in csv_reader:
            don_core.append(row[1])

    unit_list = [acc_term_L, don_core, acc_core, acc_term_R]

    return unit_list

def init_freq_lists(smiles_list):
    '''
    Creates csv files with all of the units, and sets the frequency (number of times it appears in an NFA) to 0

    Parameters
    ---------
    smiles_list: list (specific format)
        list of lists of strings
        [acceptor_terminal_left, donor_core, acceptor_core, acceptor_terminal_right]

    '''
    #freq_dir = subprocess.call('mkdir frequency_lists', shell=True)

    with open('/ihome/ghutchison/blp62/GA/running_GA/frequency_lists/freq_acc_term_updated.csv', 'w') as freq_acc_term:
        fieldnames = ['index', 'SMILES', 'frequency']
        writer = csv.DictWriter(freq_acc_term, fieldnames = fieldnames)
        writer.writeheader()
        
        for acc_term_index in range(len(smiles_list[0])):
            writer.writerow({'index': acc_term_index, 'SMILES': smiles_list[0][acc_term_index], 'frequency': 0})

    with open('/ihome/ghutchison/blp62/GA/running_GA/frequency_lists/freq_acc_core_updated.csv', 'w') as freq_acc_core:
        fieldnames = ['index', 'SMILES', 'frequency']
        writer = csv.DictWriter(freq_acc_core, fieldnames = fieldnames)
        writer.writeheader()
        
        for acc_core_index in range(len(smiles_list[2])):
            writer.writerow({'index': acc_core_index,'SMILES': smiles_list[2][acc_core_index], 'frequency': 0})

    with open('/ihome/ghutchison/blp62/GA/running_GA/frequency_lists/freq_don_core_updated.csv', 'w') as freq_don_core:
        fieldnames = ['index', 'SMILES', 'frequency']
        writer = csv.DictWriter(freq_don_core, fieldnames = fieldnames)
        writer.writeheader()
        
        for don_core_index in range(len(smiles_list[1])):
            writer.writerow({'index': don_core_index, 'SMILES': smiles_list[1][don_core_index],'frequency': 0})


def update_frequency_list(unit_type, index):
    '''
    updates the frequency csv files if the unit is used again in an NFA

    Parameters
    ---------
    unit_type: str
        options are 'acc_term', 'acc_core', or 'don_core'
    index: int
        index of the building block unit
    '''

    if unit_type == 'acc_term':
        unit_csv = "frequency_lists/freq_acc_term_updated.csv"

    elif unit_type =='acc_core':
        unit_csv = "frequency_lists/freq_acc_core_updated.csv"

    elif unit_type == 'don_core':
        unit_csv = "frequency_lists/freq_don_core_updated.csv"
    
    # reading the csv
    df = pd.read_csv(unit_csv)

    # updating frequency of unit index
    df.loc[index, 'frequency'] += 1

    # writing into the file
    df.to_csv(unit_csv, index=False)


def freq_list(gen):

    for output in glob.iglob('generations/%s/*.out' % (gen)):
            
        NFA_filename = output.split('/')[-1].split('.')[0]

        if NFA_filename not in NFAs:
            NFAs.append(NFA_filename)
            units = NFA_filename.split('_')
            update_frequency_list('acc_term', int(units[0]))
            update_frequency_list('don_core', int(units[1]))
            update_frequency_list('acc_core', int(units[2]))

smiles_list = make_unit_list()
init_freq_lists(smiles_list)

NFAs = []

for gen in range(1, 53):
    print(gen)
    freq_list(gen)

