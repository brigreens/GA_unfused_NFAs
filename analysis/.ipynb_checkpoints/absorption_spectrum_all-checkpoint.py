import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import pandas as pd
import numpy as np
import matplotlib as mpl

df = pd.read_csv('../opt_GA/GA_data_v4.csv')
df = df.sort_values(by=['PCE'], ascending=True)
df = df.reset_index(drop=True)


fig, (ax, ax2)  = plt.subplots(1, 2, figsize=(10, 4),constrained_layout=True, gridspec_kw={'width_ratios': [10, 1]})

# if performing for all GA runs
'''df_above20 = df[df.PCE > 20]
viridis = cm.get_cmap('inferno_r',len(df_above20))
newcolors = viridis(np.linspace(0, 1, len(df_above20)))'''


entry_count = len(df)
viridis = cm.get_cmap('inferno_r',50)
newcolors = viridis(np.linspace(0, 1, 50))

top_PCE = []
count = 0
for i in range(len(df)-50):
    ax.vlines(df['opt bg (nm)'][i], [0], df['oscillator strength'][i], colors = 'darkgrey')

for i in range(len(df)-50, len(df)):

    top_PCE.append(df['PCE'][i])
    ax.vlines(df['opt bg (nm)'][i], [0], df['oscillator strength'][i], colors = newcolors[count])
    count +=1

# if performing for all GA runs
'''count = 0
for i in range(len(df)-len(df_above20)):
    ax.vlines(df['opt bg (nm)'][i], [0], df['oscillator strength'][i], colors = 'darkgrey')

for i in range(len(df)-len(df_above20), len(df)):

    top_PCE.append(df['PCE'][i])
    ax.vlines(df['opt bg (nm)'][i], [0], df['oscillator strength'][i], colors = newcolors[count])
    count +=1'''

min_PCE = min(top_PCE)
max_PCE = max(top_PCE)

cmap = mpl.cm.inferno_r
norm = mpl.colors.Normalize(vmin = min_PCE, vmax = max_PCE)

#ax2 = fig.add_axes([0.1, 0.1, 0.25, 1])

#im = ax2.imshow([newcolors], orientation='vertical', cmap = viridis)
colbar = mpl.colorbar.ColorbarBase(ax2, cmap=cmap, norm=norm, orientation='vertical')
colbar.set_label('PCE', fontsize=16)
#mpl.colorbar.make_axes(ax2, aspect =5)

ax.set_xlabel('Wavelength (nm)', fontsize=16)
ax.set_ylabel('Oscillator Strength', fontsize=16)
ax.set_ylim(bottom=0)
ax.set_xlim(left=250, right=1800)

#ax2.figure.set_size_inches(2, 15)

#plt.xlabel('Wavelength (nm)')
#plt.ylabel('Oscillator Strength')

#fig.colorbar(im, ax=ax)

plt.savefig('../opt_GA/absorption_spectra_v4.png')