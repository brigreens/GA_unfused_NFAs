import csv
from openbabel import pybel
from rdkit import Chem
import glob

def RepresentsInt(s):
    try: 
        int(s)
        return True
    except ValueError:
        return False

def num_rings(SMILES):
    ring_symbol_count = 0
    frags = list(SMILES)
    print(frags)

    for x in frags:
        if RepresentsInt(x) == True:
            ring_symbol_count += 1 # adds 1 if it sees a number in the SMILES, representing a part of the ring break
        elif x == '%':
            ring_symbol_count -= 1 # the number of the ring closure is double digits, so we don't want to add that ring twice

    ring_count = ring_symbol_count / 2 # needs 2 numbers for every 1 ring break

    return ring_count

def check_mol_breaks(SMILES):
    frags = list(SMILES)
    if '.' in frags:
        return True


'''def count_rings(molecule):

    # uses Euler's rule (#Rings = #Bonds - #Atoms + #Components)
    nfrags = len(Chem.GetMolFrags(m))
    n_rings = m.GetNumBonds()-m.GetNumAtoms()+nfrags

    return n_rings'''

# not using right now
#counts = []
#counts.append(count_rings(m))
#counts.append(count_rings(opt_m))
#print(counts)

def make_SMILES_from_filename(filename, unit_list):



    filename = filename.split('/')[-1].split('.')[0]

    # splits file name by "_" symbol
    NFA_indices = filename.split('_')

    NFA_string = unit_list[0][int(NFA_indices[0])] + unit_list[1][int(NFA_indices[1])] + unit_list[2][int(NFA_indices[2])]+ unit_list[1][int(NFA_indices[1])]+ unit_list[3][int(NFA_indices[0])]

    return NFA_string


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

    with open("../running_GA/building_blocks/acc_core.csv", "r") as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        next(csv_reader)
        for row in csv_reader:
            acc_core.append(row[1])

    with open("../running_GA/building_blocks/acc_left_term.csv", "r") as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        next(csv_reader)
        for row in csv_reader:
            acc_term_L.append(row[1])

    with open("../running_GA/building_blocks/acc_right_term.csv", "r") as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        next(csv_reader)
        for row in csv_reader:
            acc_term_R.append(row[1])

    with open("../running_GA/building_blocks/don_core.csv", "r") as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        next(csv_reader)
        for row in csv_reader:
            don_core.append(row[1])

    unit_list = [acc_term_L, don_core, acc_core, acc_term_R]

    return unit_list

def xyz_to_smiles(filename):
    
    mol = next(pybel.readfile("xyz", filename))

    smi = mol.write(format="smi")

    return smi.split()[0].strip()

unit_list = make_unit_list()


for file in glob.iglob('../running_GA/opt_xyz/*.xyz'):
    filename = file.split('/')[-1].split('.')[0]
    unopt_smi = make_SMILES_from_filename(file, unit_list)
    unopt_num_rings = num_rings(unopt_smi)

    print(filename)
    opt_smi = xyz_to_smiles(file)
    opt_num_rings = num_rings(opt_smi)

    #print(opt_smi)
    if unopt_num_rings != opt_num_rings:
        print('The numbers of rings does match before and after geometry optimization')

    if check_mol_breaks(opt_smi) == True:
            print('The molecule broke into fragments')




'''for file in glob.iglob('../running_GA/opt_xyz/*.xyz'):
    filename = file.split('/')[-1].split('.')[0]
    unopt_smi = make_SMILES_from_filename(file, unit_list)
    opt_smi = xyz_to_smiles(file)
    print(filename)
    print(unopt_smi)
    try:
        canon_unopt = Chem.MolToSmiles(Chem.MolFromSmiles(unopt_smi))
        canon_opt = Chem.MolToSmiles(Chem.MolFromSmiles(opt_smi))

        if canon_unopt != canon_opt:
            
            print( filename + ' smiles did not match. Possible geometry error')

    except:
        print('an error occured with ' + filename)'''

