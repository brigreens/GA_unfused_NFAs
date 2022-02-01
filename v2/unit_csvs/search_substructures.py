from rdkit import Chem
import pandas as pd


orig_smi = 'c1c(cc2c(c1)nccn2)'
a = Chem.CanonSmiles(orig_smi)

# create list of smiles
df = pd.read_csv("don_left_term_with_homo.csv", names = ['index', 'SMILES', 'new name', 'old name', 'homo', 'lumo', 'bandgap', 'new dft filename'])
smiles = df['SMILES'].values.tolist()
smiles_don = smiles[1:]

count = 0
# search for substructures
for x in smiles_don:
    b = Chem.CanonSmiles(str(x))

    if a == b:
        print('don')
        print(orig_smi)
        print(x)
        count +=1

if count <2:
    df = pd.read_csv("acc_left_term_with_homo.csv", names = ['index', 'SMILES', 'new name', 'old name', 'homo', 'lumo', 'bandgap', 'new dft filename'])
    smiles_acc = df['SMILES'].values.tolist()
    smiles_acc = smiles_acc[1:]

    # search for substructures
    for x in smiles_acc:
        b = Chem.CanonSmiles(str(x))

        if a == b:
            print('acc')
            print(orig_smi)
            print(x)


