#path to analysis file
path = '../v1/full_analysis_data.csv'

import math
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import scipy.constants as sc
import pandas as pd

def split_filename(df):
    df2 = df.copy()

    new = df2['filename'].str.split('_', expand=True)
    df2['acc_term'] = new[0]
    # fixes issue where unit 45 and 90 are the same unit
    for x in range(len(df2['acc_term'])):
        if df2['acc_term'][x] == '45':
            df2['acc_term'][x] = '90'
    df2['don_core'] = new[1]
    df2['acc_core'] = new[2]

    return df2

# takes any dataframe with cols: acc_term, don_core, acc_core
# returns list of 3 new dataframes with cols: [acc_term, acc_term_count, acc_term_perc], [don_core, don_core_count, don_core_perc], and [acc_core, acc_core_count, acc_core_perc] 
def get_freq_counts(df):
    # split filename column
    df_split = split_filename(df)
    # make series objects of each part of filename with accompanying counts
    # note that both monomer columns are combined before creating series with value counts
    count_at = df_split['acc_term'].value_counts()
    count_dc = df_split['don_core'].value_counts()
    count_ac = df_split['acc_core'].value_counts()

    # make series objects into dataframes and rename columns appropriately
    count_at = pd.DataFrame(count_at.reset_index())
    count_at.columns = ['acc_term', 'acc_term_count']
    count_at.to_csv('acc_term_freq.csv')

    count_dc = pd.DataFrame(count_dc.reset_index())
    count_dc.columns = ['don_core', 'don_core_count']
    count_dc.to_csv('don_core_freq.csv')

    count_ac = pd.DataFrame(count_ac.reset_index())
    count_ac.columns = ['acc_core', 'acc_core_count']
    count_ac.to_csv('acc_core_freq.csv')

    # make percentage column in each dataframe   
    count_at['acc_term_perc'] = count_at['acc_term_count']/count_at['acc_term_count'].sum()*100
    count_dc['don_core_perc'] = count_dc['don_core_count']/count_dc['don_core_count'].sum()*100
    count_ac['acc_core_perc'] = count_ac['acc_core_count']/count_ac['acc_core_count'].sum()*100
    return [count_at, count_dc, count_ac]  

run_data = pd.read_csv(path, index_col=False)

at_counts, dc_counts, ac_counts = get_freq_counts(run_data)
# take top 10 only 
at_counts = at_counts[:10] 
dc_counts = dc_counts[:10] 
ac_counts = ac_counts[:10]


SMALL_SIZE = 8
MEDIUM_SIZE = 10 
BIGGER_SIZE = 12 
plt.rc('font', size=SMALL_SIZE)          # controls default text sizes 
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title 
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels 
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels 
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels 
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize 
plt.rc('figure', titlesize=MEDIUM_SIZE)  # fontsize of the figure title"


fig, axs = plt.subplots(1, 3, sharey=True, dpi=300) 

rot=45 
wt='bold' 

plt.rcParams["font.weight"] = 'regular' 

axs[0].bar('acc_term', 'acc_term_perc', data=at_counts) 
axs[0].set_title('Terminal Acceptor') 
axs[0].set_xticklabels(at_counts['acc_term'], rotation=rot) 
axs[0].set(ylim=(0, 100)) 

#axs[0].grid(alpha=0.5) 
axs[0].set_axisbelow(True) 

axs[1].bar('don_core', 'don_core_perc', data=dc_counts) 
axs[1].set_title('Core Donor') 
axs[1].set_xticklabels(dc_counts['don_core'], rotation=rot) 

#axs[1].grid(alpha=0.5) 
axs[1].set_axisbelow(True) 

for tic in axs[1].yaxis.get_major_ticks(): 
    tic.tick1On = tic.tick2On = False 

axs[2].bar('acc_core', 'acc_core_perc', data=ac_counts) 
axs[2].set_title('Core Acceptor') 
axs[2].set_xticklabels(ac_counts['acc_core'], rotation=rot) 

#axs[2].grid(alpha=0.5) 
axs[2].set_axisbelow(True) 

for tic in axs[2].yaxis.get_major_ticks(): 
    tic.tick1On = tic.tick2On = False 

top = 0.5 
bottom = 0.15 
left = 0 
right = 1 

fig.text((left+right)/2, 0.03,'Index', ha='center') 
fig.text(-0.08, (top+bottom)/2,'Frequency Percentage', va='center', rotation='vertical') 

plt.subplots_adjust(top=top, bottom=bottom, left=left, right=right, hspace=0.55, wspace=0.05) 

plt.savefig('frequency_plot.pdf', transparent=False, bbox_inches='tight') 
plt.savefig('frequency_plot.png', transparent=False, bbox_inches='tight')

#plt.show()