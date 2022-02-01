import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import pandas as pd
import numpy as np
import matplotlib as mpl
import glob
import cclib

def parse_DFT(filename): 
    """
    Parses through DFT Orca output files

    Parameters
    ----------
    filename: str
        filepath of DFT output file

    Returns
    -------
    homo_lumo: list
        homo_energy: HOMO in eV
        lumo_energy: LUMO in eV
        moment: dipole moment in units of Debye
        homo_minus1_energy: Energy of HOMO-1 in eV
        lumo_plus1_energy: Energy of LUMO+1 in eV
    """
    myfile = cclib.io.ccopen(filename)
    values = myfile.parse()
    
    homo = values.homos[0] # index of HOMO. assuming spin-restricted.
    homo_minus1 = homo -1 # index of HOMO-1 energy level
    lumo = homo + 1 # index of LUMO
    lumo_plus1 = lumo + 1 # index of LUMO + 1
    homo_energy = values.moenergies[0][homo]  #eV
    lumo_energy = values.moenergies[0][lumo]
    homo_minus1_energy = values.moenergies[0][homo_minus1]
    lumo_plus1_energy = values.moenergies[0][lumo_plus1]

    deltaHOMO = abs(homo_minus1_energy - homo_energy)
    deltaLUMO = abs(lumo_energy - lumo_plus1_energy)
    fund_bg = lumo_energy-homo_energy


    homo_lumo = [homo_minus1_energy, homo_energy, lumo_energy, lumo_plus1_energy, deltaHOMO, deltaLUMO, fund_bg]
    
    return homo_lumo


def find_new_name(oldname):
    for x in range(len(df_ac_names)):
        if oldname == df_ac_names['old name'][x]:
            newname = df_ac_names.index[x]
            outputs = [newname, 'ac']
            return outputs
    for x in range(len(df_at_names)):
        if oldname == df_at_names['old name'][x]:
            newname = df_at_names.index[x]
            outputs = [newname, 'at']
            return outputs

    for x in range(len(df_dc_names)):
        if oldname == df_dc_names['old name'][x]:
            newname = df_dc_names.index[x]
            outputs = [newname, 'dc']
            return outputs
    for x in range(len(df_dt_names)):
        if oldname == df_dt_names['old name'][x]:
            newname = df_dt_names.index[x]
            outputs = [newname, 'dt']
            return outputs

def add_freq_count(new_name, outputs, type):
    if type == 'ac':
        for x in range(len(df_ac)):
            if new_name == df_ac['acc_core'][x]:
                new_outputs = outputs.copy()
                new_outputs.append(df_ac['acc_core_count'][x])
                ac_DFT[new_name] = new_outputs
        for x in range(len(df_at)):
            if new_name == df_at['acc_term'][x]:
                new_outputs = outputs.copy()
                new_outputs.append(df_at['acc_term_count'][x])
                at_DFT[new_name] = new_outputs
                return
        outputs.append(0)
        at_DFT[new_name] = outputs
        ac_DFT[new_name] = outputs

    if type == 'at':
        for x in range(len(df_at)):
            if new_name == df_at['acc_term'][x]:
                outputs.append(df_at['acc_term_count'][x])
                at_DFT[new_name] = outputs
                return
        outputs.append(0)
        at_DFT[new_name] = outputs

    if type == 'dc':
        for x in range(len(df_dc)):
            if new_name == df_dc['don_core'][x]:
                new_outputs = outputs.copy()
                new_outputs.append(df_dc['don_core_count'][x])
                dc_DFT[new_name] = new_outputs
        for x in range(len(df_dt)):
            if new_name == df_dt['don_term'][x]:
                new_outputs = outputs.copy()
                new_outputs.append(df_dt['don_term_count'][x])
                dt_DFT[new_name] = new_outputs
                return
        outputs.append(0)
        dt_DFT[new_name] = outputs
        dc_DFT[new_name] = outputs

    if type == 'dt':
        for x in range(len(df_dt)):
            if int(new_name) == int(df_dt['don_term'][x]):
                outputs.append(df_dc['don_term_count'][x])
                dt_DFT[new_name] = outputs
                return
        outputs.append(0)
        dt_DFT[new_name] = outputs



def parse_units():
    '''
    create csv with unit energy levels. Can comment it if csv exists
    '''

    for file in glob.iglob('/ihome/ghutchison/blp62/GA/running_GA/building_blocks/optimized_B3LYP/*.out'):
        
        oldname = file.split('/')[-1].split('.',1)[0]
        # matches the new name ot the old name
        new_name_output = find_new_name(oldname)
        new_name = new_name_output[0]
        type = new_name_output[1]

        # parses the DFT file
        outputs = parse_DFT(file)

        # adds the unit frequency count to the dictionary
        add_freq_count(new_name, outputs, type)

    '''
    for file in glob.iglob('/ihome/ghutchison/blp62/GA/running_GA/more_units/new_units/*.out'):
        
        oldname = file.split('/')[-1].split('.',1)[0]
        # matches the new name ot the old name
        new_name_output = find_new_name(oldname)
        new_name = new_name_output[0]
        type = new_name_output[1]

        # parses the DFT file
        outputs = parse_DFT(file)

        # adds the unit frequency count to the dictionary
        add_freq_count(new_name, outputs, type)'''

    df_at_DFT = pd.DataFrame.from_dict(at_DFT, orient = 'index', columns = ['HOMO-1 (eV)', 'HOMO (eV)', 'LUMO (eV)', 'LUMO+1 (eV)', 'Delta HOMO', 'delta LUMO', 'fund bg', 'frequency'])
    df_at_DFT = df_at_DFT.rename_axis('Unit')
    df_at_DFT.to_csv('/ihome/ghutchison/blp62/GA/running_GA/sequences/at_freq_DFT_more_units.csv')  

    df_ac_DFT = pd.DataFrame.from_dict(ac_DFT, orient = 'index', columns = ['HOMO-1 (eV)', 'HOMO (eV)', 'LUMO (eV)', 'LUMO+1 (eV)', 'Delta HOMO', 'delta LUMO', 'fund bg', 'frequency'])
    df_ac_DFT = df_ac_DFT.rename_axis('Unit')
    df_ac_DFT.to_csv('/ihome/ghutchison/blp62/GA/running_GA/sequences/ac_freq_DFT_more_units.csv')  

    df_dc_DFT = pd.DataFrame.from_dict(dc_DFT, orient = 'index', columns = ['HOMO-1 (eV)', 'HOMO (eV)', 'LUMO (eV)', 'LUMO+1 (eV)', 'Delta HOMO', 'delta LUMO', 'fund bg', 'frequency'])
    df_dc_DFT = df_dc_DFT.rename_axis('Unit')
    df_dc_DFT.to_csv('/ihome/ghutchison/blp62/GA/running_GA/sequences/dc_freq_DFT_more_units.csv')  

    df_dt_DFT = pd.DataFrame.from_dict(dt_DFT, orient = 'index', columns = ['HOMO-1 (eV)', 'HOMO (eV)', 'LUMO (eV)', 'LUMO+1 (eV)', 'Delta HOMO', 'delta LUMO', 'fund bg', 'frequency'])
    df_dt_DFT = df_dt_DFT.rename_axis('Unit')
    df_dt_DFT.to_csv('/ihome/ghutchison/blp62/GA/running_GA/sequences/dt_freq_DFT_more_units.csv')  


at_DFT = {}
ac_DFT = {}
dc_DFT = {}
dt_DFT={}
'''
df_at_names = pd.read_csv('/ihome/ghutchison/blp62/GA/running_GA/sequence/unit_csvs/acc_left_term_fixed.csv')
df_ac_names = pd.read_csv('/ihome/ghutchison/blp62/GA/running_GA/sequence/unit_csvs/acc_core.csv')
df_dc_names = pd.read_csv('/ihome/ghutchison/blp62/GA/running_GA/sequence/unit_csvs/don_core.csv')
df_dt_names = pd.read_csv('/ihome/ghutchison/blp62/GA/running_GA/sequence/unit_csvs/don_left_term_fixed.csv')

# update with frequency cvss when GA is done
df_at = pd.read_csv('at_freq_DFT_more_units.csv')
df_at = df_at.sort_values(by=['acc_term_count'], ascending=True)
df_at = df_at.reset_index(drop=True)

df_ac = pd.read_csv('ac_freq_DFT_more_units.csv')
df_ac = df_ac.sort_values(by=['acc_core_count'], ascending=True)
df_ac = df_ac.reset_index(drop=True)

df_dc = pd.read_csv('dc_freq_DFT_more_units.csv')
df_dc = df_dc.sort_values(by=['don_core_count'], ascending=True)
df_dc = df_dc.reset_index(drop=True)

df_dt = pd.read_csv('dt_freq_DFT_more_units.csv')
df_dt = df_dt.sort_values(by=['don_term_count'], ascending=True)
df_dt = df_dtS.reset_index(drop=True)'''


df_at = pd.read_csv('/ihome/ghutchison/blp62/GA/running_GA/v1/acc_term_homo_freq.csv')
df_at = df_at.sort_values(by=['acc_term_count'], ascending=True)
df_at = df_at.reset_index(drop=True)

inferno_r_at = cm.get_cmap('inferno_r',df_at['acc_term_count'].max())
newcolors_at = inferno_r_at(np.linspace(0, 1, df_at['acc_term_count'].max()))

df_ac = pd.read_csv('/ihome/ghutchison/blp62/GA/running_GA/v1/acc_core_homo_freq.csv')
df_ac = df_ac.sort_values(by=['acc_core_count'], ascending=True)
df_ac = df_ac.reset_index(drop=True)

inferno_r_ac = cm.get_cmap('inferno_r',df_ac['acc_core_count'].max())
newcolors_ac = inferno_r_ac(np.linspace(0, 1, df_ac['acc_core_count'].max()))

df_dc = pd.read_csv('/ihome/ghutchison/blp62/GA/running_GA/v1/don_core_homo_freq.csv')
df_dc = df_dc.sort_values(by=['don_core_count'], ascending=True)
df_dc = df_dc.reset_index(drop=True)

inferno_r_dc = cm.get_cmap('inferno_r',df_dc['don_core_count'].max())
newcolors_dc = inferno_r_dc(np.linspace(0, 1, df_dc['don_core_count'].max()))

'''
df_dt = pd.read_csv('dt_freq_DFT_more_units.csv')
df_dt = df_dt.sort_values(by=['frequency'], ascending=True)
df_dt = df_dtS.reset_index(drop=True)'''

#parse_units()   

#entry_count = len(df)
#viridis = cm.get_cmap('viridis',entry_count)
#newcolors = viridis(np.linspace(0, 1, entry_count))

def plot_homolumo(df, column_name, colors, ax, cb_ax, title):
    entry_count = len(df)
    top_freq = []
    count = 0

    print(df)
    for i in range(entry_count):
        ax.scatter(x=df['bandgap'][i], y=df['lumo'][i], c=[colors[df[column_name][i]-1]])
        ax.scatter(x=df['bandgap'][i], y=df['homo'][i], c=[colors[df[column_name][i]-1]])
        count +=1
        top_freq.append(df[column_name][i])

    ax.set_xlabel('Fundamental Bandgap (eV)', fontsize=16)
    ax.set_ylabel('Energy (eV)', fontsize=16)
    ax.set_xlim(1, 9)
    ax.set_ylim(-10, 3)
    ax.title.set_text(title)
    ax.title.set_size(20)
    ax.text(x=2, y=-9, s= 'HOMO', fontsize=14)
    ax.text(x=2, y=-1, s= 'LUMO', fontsize=14)

    min_freq = min(top_freq)
    max_freq = max(top_freq)

    cmap = mpl.cm.inferno_r
    norm = mpl.colors.Normalize(vmin = min_freq, vmax = max_freq)

    colbar = mpl.colorbar.ColorbarBase(cb_ax, cmap=cmap, norm=norm, orientation='vertical')
    colbar.set_label('Frequency', fontsize=16)


fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6))  = plt.subplots(3, 2, figsize=(10, 10),constrained_layout=True, gridspec_kw={'width_ratios': [10, 1]})


plot_homolumo(df_at, 'acc_term_count', newcolors_at, ax1, ax2, 'Acceptor Terminal Units')
plot_homolumo(df_ac, 'acc_core_count', newcolors_ac, ax3, ax4, 'Acceptor Core Units')
plot_homolumo(df_dc, 'don_core_count', newcolors_dc, ax5, ax6, 'Donor Core Units')


#ax.axes.xaxis.set_visible(False)

plt.savefig('homo_lumo_units_v1.png')
plt.savefig('homo_lumo_units_v1.pdf')