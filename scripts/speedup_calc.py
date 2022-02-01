# calculate speedup of GA vs brute force

import glob

# number of generations required to get to highest PCE
convg_gen = 84

# total number of possible molecules in search space
search_space = 2.05 * 10**15

# go through generations to count how many unique molecules were analyzed
unique_mol = []
for gen in range(1, convg_gen+1):
    for file in glob.iglob('/ihome/ghutchison/blp62/GA/running_GA/v5/generations/%s/*.out' % (str(gen))):
        mol = file.split('/')[-1].split('.')[0]
        if mol not in unique_mol:
            unique_mol.append(mol)

speedup = search_space / len(unique_mol)

print('Total number of molecules analyzed is ' + str(len(unique_mol)))

print('Total speedup over brute force is ' + str(speedup))