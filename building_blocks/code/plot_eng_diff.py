
import pandas as pd
import matplotlib.cm as cm
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import seaborn as sns


# open file with energy diffence compared to thiophene
with open('eng_diff.txt', 'r') as file:
    eng_diff_str = file.readlines()

# convert string txt to floats
eng_diff_str = eng_diff_str[0].split("[")[-1].split(']')[0].split(',')
energy_difference  = [float(i)*-1 for i in eng_diff_str]
# sort in increasing order
energy_difference.sort()

N = len(energy_difference)
x = range(len(energy_difference))

data = pd.DataFrame({'Eng_diff': energy_difference, 'Unit_Index': x})

# sets the diverging color palette. 245 = blue, 0 = red
colors_pal = sns.diverging_palette(245, 0, n= N, sep = 10)

fig, ax = plt.subplots()

sns.barplot(x='Unit_Index', y='Eng_diff', data = data, palette = colors_pal, ax = ax)
#ax.set_title('Identification of Electron Donor/Acceptor Units', fontsize = 16)
ax.set_xlabel('Building Block Index', fontsize = 14)
ax.set_ylabel('Building Block HOMO - Thiophene HOMO (eV)', fontsize = 12)
# remove x axis ticks and tick labels
ax.tick_params(axis = 'x', which = 'both', labelbottom = False, bottom=False)

# inverse y-axis so it goes from positive to negative


ax.legend(['Electron Withdrawing', 'Electron Donating'], loc = 'upper left', fontsize = 14)
leg = ax.get_legend()
leg.legendHandles[1].set_color('palevioletred')
leg.legendHandles[0].set_color('steelblue')

plt.savefig('distribution_of_HOMO_updated.png')
plt.savefig('distribution_of_HOMO_updated.pdf')

