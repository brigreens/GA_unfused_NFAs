import csv
import glob
import pandas as pd
from rdkit import Chem
from openbabel import pybel
from rdkit.Chem import Draw

def make3D(mol):
    '''
    Makes the mol object from SMILES 3D

    Parameters
    ---------
    mol: object
        pybel molecule object
    '''
    # make mol object 3D and add hydrogens
    pybel._builder.Build(mol.OBMol)
    mol.addh()

    ff = pybel._forcefields["mmff94"]
    success = ff.Setup(mol.OBMol)
    if not success:
        ff = pybel._forcefields["uff"]
        success = ff.Setup(mol.OBMol)
        if not success:
            sys.exit("Cannot set up forcefield")

    ff.ConjugateGradients(100, 1.0e-3)
    ff.WeightedRotorSearch(100, 25)
    ff.ConjugateGradients(250, 1.0e-4)

    ff.GetCoordinates(mol.OBMol)

def RepresentsInt(s):
    '''
    Checks if character is an integer

    Parameters
    ----------
    s: str
        character in string

    Returns
    -------
    True if it is an integer, False if it is not
    '''
    try: 
        int(s)
        return True
    except ValueError:
        return False

def num_rings(SMILES):
    '''
    Counts the number of rings in SMILES

    Parameters
    ----------
    SMILES: str
        SMILES string of the molecule

    Returns
    --------
    ring_count: int
        number of rings in the molecule

    '''

    ring_symbol_count = 0
    frags = list(SMILES)

    for x in frags:
        if RepresentsInt(x) == True:
            ring_symbol_count += 1 # adds 1 if it sees a number in the SMILES, representing a part of the ring break
        elif x == '%':
            ring_symbol_count -= 1 # the number of the ring closure is double digits, so we don't want to add that ring twice

    ring_count = ring_symbol_count / 2 # needs 2 numbers for every 1 ring break

    return ring_count


# create list of current canonical SMILES so we prevent duplication

df_orig = pd.read_csv("original_unit_SMILES.csv", names = ['SMILES'])
orig_smiles = df_orig['SMILES'].values.tolist()

unique_smi = []
mols = []
ind = []

temp_ind = 0
for x in orig_smiles[1:]:
    if num_rings(x) == 0:
        mol = Chem.MolFromSmiles(x)
        mols.append(mol)
        ind.append(str(temp_ind))
    temp_ind +=1

    canonical_smiles = Chem.MolToSmiles(Chem.MolFromSmiles(x))
    if canonical_smiles not in unique_smi:
        unique_smi.append(canonical_smiles)
    else:
        print('duplicate smile is: ' + str(x))
        index = unique_smi.index(canonical_smiles)
        print('Original smile is at index ' + str(index) + ' and is: ' + unique_smi[index])

    




# create new xyz files for new monomers
with open('1235MonomerList.txt') as new_mons:

    count = 0
    line_count = 1
    for line in new_mons:

        if num_rings(x) == 0:
            mol = Chem.MolFromSmiles(x)
            mols.append(mol)
            ind.append(str(temp_ind))
        temp_ind +=1

        canonical_smiles = Chem.MolToSmiles(Chem.MolFromSmiles(line))
        if canonical_smiles not in unique_smi:
            unique_smi.append(canonical_smiles)
            #mol = pybel.readstring('smi', line)
            #make3D(mol)
            #mol.write('xyz', '/ihome/ghutchison/blp62/GA/running_GA/more_units/new_units/%s.xyz' % (count), overwrite=False)
            #append_list_as_row('/ihome/ghutchison/blp62/GA/running_GA/more_units/unit_csvs/matching_filename_numbers_2.csv', [line_count, count, line])

            count +=1

        else:
            #append_list_as_row('/ihome/ghutchison/blp62/GA/running_GA/more_units/unit_csvs/matching_filename_numbers_2.csv', [line_count, 'N/A', line])
            #print('duplicate ' + str(count))
            #print(line)
            
            print('duplicate smile is: ' + str(x))
            index = unique_smi.index(canonical_smiles)
            print('Original smile is at index ' + str(index) + ' and is: ' + unique_smi[index])

            #print('Original smiles is ' + str(unique_smi[index] + ' located at index ' + str(index)))

            count +=1 

        line_count +=1
            
def save_svg(svg, svg_file, dpi=300):
    png_file = svg_file.replace('.svg', '.png')
    with open(svg_file, 'w') as afile:
        afile.write(svg)
    a_str = svg.encode('utf-8')
    return
          
#img = Draw.MolsToGridImage(mols, molsPerRow=4, legends = ind, subImgSize=(250, 250), returnPNG=False, useSVG=True)
#save_svg(img, 'drawing_molecules_0rings.svg', dpi=200)

def append_list_as_row(file_name, list_of_elem):
    # Open file in append mode
    with open(file_name, 'a+', newline='') as write_obj: 
        # Create a writer object from csv module
        csv_writer = csv.writer(write_obj, quoting = csv.QUOTE_NONE, escapechar=',')
        # Add contents of list as last row in the csv file
        csv_writer.writerow(list_of_elem)