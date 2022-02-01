import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import matplotlib.colors as colors
import glob
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import matplotlib as mpl


df_acc = pd.read_csv('/ihome/ghutchison/blp62/GA/running_GA/v5/unit_csvs/acc_left_term_v5.csv')
df_don = pd.read_csv('/ihome/ghutchison/blp62/GA/running_GA/v5/unit_csvs/don_left_term_v5.csv')

def even_or_odd(num):
    '''
    determines if number is even or odd

    Parameters
    --------
    num: int
        number to check if even or odd

    Returns
    -------
    'even' is number is even, 'odd' if number is odd
    '''
    if int(num) % 2 == 0:
        return 'even'
    elif int(num) == 0:
        return 'even'
    else:
        return 'odd'

HOMO_acc = df_acc['homo']
HOMO_don = df_don['homo']

def find_HOMOs(filename):

    print(filename)
    sequence = filename.split('_')[-1]
    seq = list(sequence)
    units = filename.split('_')[:-1]
    HOMOs = []

    for x in range(5):
        unit = units[x]

        if even_or_odd(seq[x]) == 'even':
            # donor
            HOMOs.append(HOMO_don[int(unit)])
        else:
            # acceptor
            HOMOs.append(HOMO_acc[int(unit)])

    return HOMOs


class MidpointNormalize(colors.Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y), np.isnan(value))

def plot_barcode(ax, values):

    rows = 1
    cols = 5

    x = np.arange(cols + 1)
    y = np.arange(rows + 1)

    Z = np.array(values).reshape(rows, cols)

    #ax.pcolormesh(x, y, Z, shading='flat', cmap='bwr_r', norm=MidpointNormalize(midpoint=-6.5467, vmin=Z.min(), vmax = Z.max()))
    ax.pcolormesh(x, y, Z, shading='flat', cmap='bwr_r', norm=MidpointNormalize(midpoint=-6.5467, vmin=-8.89652, vmax = -3.73493))

    for val in range(5):
        ax.text(x = val + 0.5, y = 0.5, s=str(round(values[val], 2)), horizontalalignment='center', verticalalignment='center', fontsize = 16)

    ax.axes.get_xaxis().set_visible(False)
    ax.axes.get_yaxis().set_visible(False)



fig = plt.figure(figsize=(20, 12), constrained_layout=True)
gs = fig.add_gridspec(ncols = 17, nrows = 9)
ax_0 = fig.add_subplot(gs[0, :5])
ax_1= fig.add_subplot(gs[1, :5])
ax_2 = fig.add_subplot(gs[2, :5])
ax_3 = fig.add_subplot(gs[3, :5])
ax_4 = fig.add_subplot(gs[4, :5])
ax_5 = fig.add_subplot(gs[5, :5])
ax_6 = fig.add_subplot(gs[6, :5])
ax_7 = fig.add_subplot(gs[7, :5])
ax_8 = fig.add_subplot(gs[8, :5])

ax_9 = fig.add_subplot(gs[0, 6:11])
ax_10 = fig.add_subplot(gs[1, 6:11])
ax_11= fig.add_subplot(gs[2, 6:11])
ax_12= fig.add_subplot(gs[3, 6:11])
ax_13 = fig.add_subplot(gs[4, 6:11])
ax_14= fig.add_subplot(gs[5, 6:11])
ax_15= fig.add_subplot(gs[6, 6:11])
ax_16= fig.add_subplot(gs[7, 6:11])
ax_17= fig.add_subplot(gs[8, 6:11])

ax_18= fig.add_subplot(gs[0, 12:])
ax_19= fig.add_subplot(gs[1, 12:])
ax_20= fig.add_subplot(gs[2, 12:])
ax_21= fig.add_subplot(gs[3, 12:])
ax_22= fig.add_subplot(gs[4, 12:])
ax_23= fig.add_subplot(gs[5, 12:])
ax_24= fig.add_subplot(gs[6, 12:])

ax_25 = fig.add_subplot(gs[7, 12:])



axes = [ax_0, ax_1, ax_2, ax_3, ax_4, ax_5, ax_6, ax_7, ax_8, ax_9, ax_10, ax_11, ax_12, ax_13, ax_14, ax_15, ax_16, ax_17, ax_18, ax_19, ax_20, ax_21, ax_22, ax_23, ax_24]

top_NFAs = pd.read_csv('/ihome/ghutchison/blp62/GA/running_GA/v5/top25.csv')
top_NFAs_filenames = top_NFAs['filename'][:25]
top_PCEs = top_NFAs['PCE'][:25]

'''
count = 0
for file in glob.iglob('/ihome/ghutchison/blp62/GA/running_GA/sequence/top10/*.out'):
    values = find_HOMOs(file)
    plot_barcode(axes[count], values)

    filename = file.split('/')[-1].split('.',1)[0]

    axes[count].set_title(filename, fontsize = 8)

    count +=1'''

# uncomment this section to make top 25 barcode figure
'''for i in range(len(top_NFAs_filenames)):
    values = find_HOMOs(top_NFAs_filenames[i])
    plot_barcode(axes[i], values)

    PCE = str(round(top_PCEs[i], 2))
    filename = top_NFAs_filenames[i].split('/')[-1].split('.',1)[0]
    title = filename + ', PCE = ' + PCE + '%'

    axes[i].set_title(title, fontsize = 16)

plt.savefig('barcode_v5.png')
plt.savefig('barcode_v5.pdf')'''


'''
Use this to plot barcodes of specific files
'''
fig = plt.figure(figsize=(8, 12), constrained_layout=True)

filenames = ['281_34_580_22_281_13051', '259_34_580_22_259_13051', '112_34_580_22_112_13051', '6_34_580_22_6_13051', '129_34_580_22_129_13051', '66_34_580_22_66_13051', '7_34_580_22_7_13051','154_34_580_22_154_13051', '192_34_580_22_192_13051', '245_34_580_22_245_13051' ]
PCEs = [20.2362, 15.38477, 15.09746, 14.92827, 14.86825, 14.84831, 14.6302, 14.4139, 14.10718, 13.92275]
num_files = len(filenames)
print(num_files)
print(len(PCEs))

for i in range(1, num_files+1):
    ax = fig.add_subplot(num_files, 1, i)
    
    values = find_HOMOs(filenames[i-1])
    plot_barcode(ax, values)

    PCE = str(round(PCEs[i-1], 2))
    title = filenames[i-1] + ', PCE = ' + PCE + '%'

    ax.set_title(title, fontsize = 16)

fig.subplots_adjust(hspace=0.6)

plt.savefig('barcode_v5_diff_acc_term.png')
plt.savefig('barcode_v5_diff_acc_term.pdf')