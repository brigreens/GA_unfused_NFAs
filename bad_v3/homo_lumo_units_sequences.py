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


def find_new_name(oldname, type):

    if type == 'acc':
        '''
        for x in range(len(df_acc_core_names)):
            if oldname == df_acc_core_names['old name'][x]:
                newname = df_acc_core_names.index[x]
                outputs = [newname, 'ac']
                return outputs
        '''

        for x in range(len(df_acc_names)):
            if oldname == df_acc_names['old name'][x]:
                newname = df_acc_names.index[x]
                outputs = [newname, 'at']
                return outputs

    if type == 'don':
        for x in range(len(df_don_names)):
            if oldname == df_don_names['old name'][x]:
                newname = df_don_names.index[x]
                outputs = [newname, 'dt']
                return outputs

def add_freq_count(new_name, outputs, type):

    '''if type == 'ac':
        check = False
        for x in range(len(df_ac)):
            if new_name == df_ac['acc_core'][x]:
                new_outputs = outputs.copy()
                new_outputs.append(df_ac['acc_core_count'][x])
                ac_DFT[new_name] = new_outputs
                check=True
        # did not find it in the acc cores
        if check == False:
            outputs.append(0)
            ac_DFT[new_name] = outputs
            at_DFT[new_name] = outputs
            return

        for x in range(len(df_at)):
            if new_name == df_at['acc_term'][x]:
                new_outputs = outputs.copy()
                new_outputs.append(df_at['acc_term_count'][x])
                at_DFT[new_name] = new_outputs
                return
        outputs.append(0)
        at_DFT[new_name] = outputs
        ac_DFT[new_name] = outputs'''


    if type == 'at':
        for x in range(len(df_at)):
            if int(new_name) == int(df_at['acc_term'][x]):
                outputs.append(df_at['acc_term_count'][x])
                at_DFT[new_name] = outputs
                return
        outputs.append(0)
        at_DFT[new_name] = outputs
                
        
    if type == 'dt':
        for x in range(len(df_dc)):

            if int(new_name) == int(df_dc['don_core'][x]):
                outputs.append(df_dc['don_core_count'][x])
                dt_DFT[new_name] = outputs
                return
            
        outputs.append(0)
        dt_DFT[new_name] = outputs

def parse_units():
    '''
    create csv with unit energy levels. Can comment it if csv exists
    '''

    for file in glob.iglob('/ihome/ghutchison/blp62/GA/running_GA/building_blocks/optimized_B3LYP/*.out'):
        
        # parses the DFT file
        outputs = parse_DFT(file)

        oldname = file.split('/')[-1].split('.',1)[0]
        
        # matches the new name to the old name
        if outputs[1] < -6.5467:
            new_name_output = find_new_name(oldname, 'acc')
        else:  
            new_name_output = find_new_name(oldname, 'don')

        new_name = new_name_output[0]
        type = new_name_output[1]

        # adds the unit frequency count to the dictionary
        add_freq_count(new_name, outputs, type)

    
    for file in glob.iglob('/ihome/ghutchison/blp62/GA/running_GA/more_units/new_units/*.out'):
        
        # parses the DFT file
        outputs = parse_DFT(file)

        oldname = file.split('/')[-1].split('.',1)[0]
        
        # matches the new name to the old name
        if outputs[1] < -6.5467:
            new_name_output = find_new_name(oldname, 'acc')
        else:  
            new_name_output = find_new_name(oldname, 'don')

        new_name = new_name_output[0]
        type = new_name_output[1]

        # adds the unit frequency count to the dictionary
        add_freq_count(new_name, outputs, type)
    

    df_acc_DFT = pd.DataFrame.from_dict(at_DFT, orient = 'index', columns = ['HOMO-1 (eV)', 'HOMO (eV)', 'LUMO (eV)', 'LUMO+1 (eV)', 'Delta HOMO', 'delta LUMO', 'fund bg', 'frequency'])
    df_acc_DFT = df_acc_DFT.rename_axis('Unit')
    df_acc_DFT.to_csv('/ihome/ghutchison/blp62/GA/running_GA/more_units/acc_freq_DFT_more_units.csv')  

    df_don_DFT = pd.DataFrame.from_dict(dt_DFT, orient = 'index', columns = ['HOMO-1 (eV)', 'HOMO (eV)', 'LUMO (eV)', 'LUMO+1 (eV)', 'Delta HOMO', 'delta LUMO', 'fund bg', 'frequency'])
    df_don_DFT = df_don_DFT.rename_axis('Unit')
    df_don_DFT.to_csv('/ihome/ghutchison/blp62/GA/running_GA/more_units/don_freq_DFT_more_units.csv')  


at_DFT = {}
dt_DFT={}
ac_DFT={}

'''df_acc_names = pd.read_csv('/ihome/ghutchison/blp62/GA/running_GA/sequence/unit_csvs/acc_left_term_fixed.csv')
df_don_names = pd.read_csv('/ihome/ghutchison/blp62/GA/running_GA/sequence/unit_csvs/don_left_term_fixed.csv')
'''
df_acc_names = pd.read_csv('/ihome/ghutchison/blp62/GA/running_GA/more_units/unit_csvs/acc_left_term.csv')
df_don_names = pd.read_csv('/ihome/ghutchison/blp62/GA/running_GA/more_units/unit_csvs/don_left_term.csv')
df_acc_core_names = pd.read_csv('/ihome/ghutchison/blp62/GA/running_GA/more_units/unit_csvs/acc_core.csv')


'''df_acc = pd.read_csv('acc_unit_freq_seq.csv')
df_acc = df_acc.sort_values(by=['acc_unit_count'], ascending=True)
df_acc = df_acc.reset_index(drop=True)

df_don = pd.read_csv('don_unit_freq_seq.csv')
df_don = df_don.sort_values(by=['don_unit_count'], ascending=True)
df_don = df_don.reset_index(drop=True)'''

df_ac = pd.read_csv('../analysis/acc_core_freq_more_units.csv')
df_ac = df_ac.sort_values(by=['acc_core_count'], ascending=True)
df_ac = df_ac.reset_index(drop=True)

df_at = pd.read_csv('../analysis/acc_term_freq_more_units.csv')
df_at = df_at.sort_values(by=['acc_term_count'], ascending=True)
df_at = df_at.reset_index(drop=True)

df_dc = pd.read_csv('../analysis/don_core_freq_more_units.csv')
df_dc = df_dc.sort_values(by=['don_core_count'], ascending=True)
df_dc = df_dc.reset_index(drop=True)

print(df_ac)
print(len(df_ac))

print(df_at)
print(len(df_at))

print(df_dc)
print(len(df_dc))

'''df_acc = pd.read_csv('at_freq_DFT_more_units.csv')
df_acc = df_acc.sort_values(by=['frequency'], ascending=True)
df_acc = df_acc.reset_index(drop=True)

viridis_at = cm.get_cmap('autumn_r',668)
newcolors_at = viridis_at(np.linspace(0, 1, 668))

df_ac = pd.read_csv('ac_freq_DFT_more_units.csv')
df_ac = df_ac.sort_values(by=['frequency'], ascending=True)
df_ac = df_ac.reset_index(drop=True)

viridis_ac = cm.get_cmap('autumn_r',1814)
newcolors_ac = viridis_ac(np.linspace(0, 1, 1814))

df_dc = pd.read_csv('dc_freq_DFT_more_units.csv')
df_dc = df_dc.sort_values(by=['frequency'], ascending=True)
df_dc = df_dc.reset_index(drop=True)

viridis_dc = cm.get_cmap('autumn_r',588)
newcolors_dc = viridis_dc(np.linspace(0, 1, 588))


df_don = pd.read_csv('dt_freq_DFT_more_units.csv')
df_don = df_don.sort_values(by=['frequency'], ascending=True)
df_don = df_don.reset_index(drop=True)'''

parse_units()   

#entry_count = len(df)
#viridis = cm.get_cmap('viridis',entry_count)
#newcolors = viridis(np.linspace(0, 1, entry_count))

def plot_homolumo(df, column_name, colors, ax, cb_ax, title):
    entry_count = len(df)
    top_freq = []
    count = 0


    for i in range(entry_count):
        ax.scatter(x=df['fund bg'][i], y=df['LUMO (eV)'][i], c=[colors[df['frequency'][i]-1]])
        ax.scatter(x=df['fund bg'][i], y=df['HOMO (eV)'][i], c=[colors[df['frequency'][i]-1]])
        count +=1
        top_freq.append(df[column_name][i])

    ax.set_xlabel('Fundamental Bandgap (eV)', fontsize=16)
    ax.set_ylabel('Energy (eV)', fontsize=16)
    ax.set_xlim(1, 9)
    ax.set_ylim(-10, 3)
    ax.title.set_text(title)
    ax.title.set_size(20)
    

    min_freq = min(top_freq)
    max_freq = max(top_freq)

    cmap = mpl.cm.autumn_r
    norm = mpl.colors.Normalize(vmin = min_freq, vmax = max_freq)

    colbar = mpl.colorbar.ColorbarBase(cb_ax, cmap=cmap, norm=norm, orientation='vertical')
    colbar.set_label('Frequency', fontsize=16)
'''

fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6))  = plt.subplots(3, 2, figsize=(10, 10),constrained_layout=True, gridspec_kw={'width_ratios': [10, 1]})


plot_homolumo(df_acc, 'frequency', newcolors_at, ax1, ax2, 'Acceptor Terminal Units')
plot_homolumo(df_ac, 'frequency', newcolors_ac, ax3, ax4, 'Acceptor Core Units')
plot_homolumo(df_dc, 'frequency', newcolors_dc, ax5, ax6, 'Donor Core Units')


#ax.axes.xaxis.set_visible(False)

plt.savefig('homo_lumo_building_blocks_more_units.png')'''