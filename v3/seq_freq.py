
#path to analysis file
path = 'full_analysis_data.csv'

import math
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import scipy.constants as sc
import pandas as pd

df_acc = pd.read_csv('/ihome/ghutchison/blp62/GA/running_GA/v5/unit_csvs/acc_left_term_v5.csv')
df_don = pd.read_csv('/ihome/ghutchison/blp62/GA/running_GA/v5/unit_csvs/don_left_term_v5.csv')


def sym(df):

    syms = []

    for x in df['sequence']:

        # finding the mid, start and last index of the string
        mid = (len(x)-1)/2
        start = 0
        last = len(x)-1
        flag = 0
    
        # A loop until the mid of the string
        while(start<mid):
    
            # comparing letters from right from the letters from left
            if (x[start]== x[last]):
                
                start += 1
                last -= 1
                
            else:
                flag = 1
                break
                
        # Checking the flag variable to check if the string is palindrome or not
        if flag == 0:
            syms.append('Symmetrical')
        else:
            syms.append('Asymmetrical')

    return syms

def unit_freq(df):
    df3 = df.copy()

    new = df3['filename'].str.split('_', expand=True)
    print(new)
    df3['sequence'] = new[new.columns[-1]]
    
    acceptors = []
    donors = []
    count = 0
    for x in range(len(df3['sequence'])):
        seq_list = list(df3['sequence'][x])
        for i in range(len(seq_list)):
            if even_or_odd(seq_list[i]) == 'odd':
                #print(new[i][count])
                acceptors.append(new[i][count])

            else:
                donors.append(new[i][count])
        count +=1

    df_acc = pd.DataFrame(acceptors, columns=['acc units'])
    df_don = pd.DataFrame(donors, columns=['don units'])

    count_acc = df_acc['acc units'].value_counts()
    count_don = df_don['don units'].value_counts()

    # make series objects into dataframes and rename columns appropriately
    count_acc = pd.DataFrame(count_acc.reset_index())
    count_acc.columns = ['acc_unit', 'acc_unit_count']
    count_acc.to_csv('acc_unit_freq_seq.csv')

    count_don = pd.DataFrame(count_don.reset_index())
    count_don.columns = ['don_unit', 'don_unit_count']
    count_don.to_csv('don_unit_freq_seq.csv')

    # make percentage column in each dataframe   
    count_acc['acc_unit_perc'] = count_acc['acc_unit_count']/count_acc['acc_unit_count'].sum()*100
    count_don['don_unit_perc'] = count_don['don_unit_count']/count_don['don_unit_count'].sum()*100

    acc_counts = count_acc[:10] 
    don_counts = count_don[:10] 

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


    fig, axs = plt.subplots(1, 2, sharey=True, dpi=300) 

    rot=45 
    wt='bold' 

    plt.rcParams["font.weight"] = 'regular' 

    axs[0].bar('acc_unit', 'acc_unit_perc', data=acc_counts) 
    axs[0].set_title('Acceptor Units') 
    axs[0].set_xticklabels(acc_counts['acc_unit'], rotation=rot) 
    axs[0].set(ylim=(0, 100)) 
    #axs[0].grid(alpha=0.5) 
    axs[0].set_axisbelow(True) 

    axs[1].bar('don_unit', 'don_unit_perc', data=don_counts) 
    axs[1].set_title('Donor Units') 
    axs[1].set_xticklabels(don_counts['don_unit'], rotation=rot) 
    #axs[1].grid(alpha=0.5) 
    axs[1].set_axisbelow(True) 
    for tic in axs[1].yaxis.get_major_ticks(): 
        tic.tick1On = tic.tick2On = False 


    top = 0.5 
    bottom = 0.15 
    left = 0 
    right = 1 

    fig.text((left+right)/2, 0.03,'Index', ha='center') 
    fig.text(-0.08, (top+bottom)/2,'Frequency Percentage', va='center', rotation='vertical') 

    plt.subplots_adjust(top=top, bottom=bottom, left=left, right=right, hspace=0.55, wspace=0.05) 

    plt.savefig('unit_frequency_plot.pdf', transparent=False, bbox_inches='tight') 
    plt.savefig('unit_frequency_plot.png', transparent=False, bbox_inches='tight')

def find_homo(unit_index, row, new_df, unit_AD):
        
        # index label of the first unit
        unit = int(row[new_df.columns[unit_index]])

        # if unit1_AD is even, it is a donor so look at donor csv
        if even_or_odd(unit_AD) == 'even':
            homo = df_don['homo'][unit]
        # if unit1_AD is odd, it is an acceptor so look at acc csv
        else:
            homo = df_acc['homo'][unit]

        return homo

def sequence(new_df):
    '''
    Finds the actual sequence (fixes error from mismatched DFT files) by comparing the HOMO directly to thiophene
    '''
    sequences = []
    # iterate through dataframe to update the sequences
    for index, row in new_df.iterrows():
        # initialize empty new sequence
        new_seq = ''
        # intialize acc or don count
        don_count = 0
        acc_count = 1

        seq_tracker = []
        # string of the original sequence, so we know how the units were classified 
        orig_seq = str(row[new_df.columns[-1]])

        for x in range(5):
            # find number (0 or 1) of the first unit for sequence
            unit_AD = orig_seq[x]

            homo = find_homo(x, row, new_df, unit_AD)

            if unit_AD in seq_tracker:
                temp_AD_index = int(seq_tracker.index(unit_AD))
                new_seq += str(new_seq[temp_AD_index])

            # acceptor
            elif homo < -6.5467:
                new_seq += str(acc_count)
                acc_count +=2
            # donor
            else:
                new_seq += str(don_count)
                don_count +=2

            seq_tracker.append(unit_AD)

        sequences.append(new_seq)
        #new_df.loc[index, 'sequence'] = new_seq

    new_df['updated_sequence'] = sequences

    return new_df

    

def split_filename(df):
    '''
    Takes the full_analysis_data and adds a column for sequence and symmetry
    Parameters
    --------
    df: dataframe

    Returns
    -------
    df2: dataframe
        new dataframe with sequence and symmetry added
    '''

    df2 = df.copy()

    new = df2['filename'].str.split('_', expand=True)

    # finds the updated sequence and adds it to the new df
    sequence(new)
    #seqs = list(new)
    df2['sequence'] = new[new.columns[-1]]
    
    df2['sym'] = sym(df2)

    return df2

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

def useAandD(count_seq):
    seq_list = count_seq['sequence']
    new_seq_strs = []

    for x in seq_list:
        seq_str = list(x)

        new_seq = []
        for i in seq_str:
            if i == '0':
                new_seq.append('D')
            elif i == '1':
                new_seq.append('A')
            elif i == '2':
                new_seq.append("D'")
            elif i == '3':
                new_seq.append("A'")
            elif i == '4':
                new_seq.append('D"')
            elif i == '5':
                new_seq.append('A"')
            elif i == '6':
                new_seq.append('D"\'')
            elif i == '7':
                new_seq.append('A"\'')
            elif i == '8':
                new_seq.append('D""')
            elif i == '9':
                new_seq.append('A""')

        new_seq_str = new_seq[0] + '-' + new_seq[1] + '-' + new_seq[2] + '-' + new_seq[3] + '-' + new_seq[4] 
        new_seq_strs.append(new_seq_str)
        #print(new_seq_str)

    return new_seq_strs


def get_freq_counts(df):
    # split filename column
    print(df)
    df_split = split_filename(df)
    print(df_split)
    # make series objects of each part of filename with accompanying counts
    # note that both monomer columns are combined before creating series with value counts
    count_seq = df_split['sequence'].value_counts()

    count_sym = df_split['sym'].value_counts()

    # make series objects into dataframes and rename columns appropriately
    count_seq = pd.DataFrame(count_seq.reset_index())
    count_seq.columns = ['sequence', 'sequence_count']

    count_sym = pd.DataFrame(count_sym.reset_index())
    count_sym.columns = ['Symmetry', 'sym_count']

    # make percentage column in each dataframe   
    count_seq['sequence_count_perc'] = count_seq['sequence_count']/count_seq['sequence_count'].sum()*100

    count_seq['sequence_AD'] = useAandD(count_seq)

    count_sym['sym_count_perc'] = count_sym['sym_count']/count_sym['sym_count'].sum()*100

    
    return [count_seq, count_sym], df_split

run_data = pd.read_csv(path, index_col=False)

counts, df_updated = get_freq_counts(run_data)
#print(counts)
print(df_updated)

#this creates plot for frequency of units. This will be incorrect and its not necessary
#unit_freq(run_data)


seq_count = counts[0]
sym_count = counts[1]

# take top 10 only 
at_counts = seq_count[:10]



SMALL_SIZE = 8
MEDIUM_SIZE = 10 
BIGGER_SIZE = 12 
plt.rc('font', size=BIGGER_SIZE)          # controls default text sizes 
plt.rc('axes', titlesize=BIGGER_SIZE)     # fontsize of the axes title 
plt.rc('axes', labelsize=BIGGER_SIZE)    # fontsize of the x and y labels 
plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels 
plt.rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels 
plt.rc('legend', fontsize=MEDIUM_SIZE)    # legend fontsize 
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title"


fig, ax = plt.subplots(dpi=300) 

rot=45 
wt='bold' 

plt.rcParams["font.weight"] = 'regular' 

ax.bar('sequence_AD', 'sequence_count_perc', data=at_counts) 
#ax.set_title('Sequence') 
ax.set_xticklabels(seq_count['sequence_AD'], rotation=rot) 
ax.set(ylim=(0, 100)) 

ax.grid(axis = 'y', alpha=0.5) 
ax.set_axisbelow(True) 

ax.set_ylabel('Frequency Percentage')
#ax.set_xlabel('Sequence')

plt.savefig('../v5/seq_frequency_v5.pdf', transparent=False, bbox_inches='tight') 
plt.savefig('../v5/seq_frequency_v5.png', transparent=False, bbox_inches='tight')

fig, ax2 = plt.subplots(dpi=300) 

ax2.bar('Symmetry', 'sym_count_perc', data=sym_count) 
#ax2.set_title('Symmetry') 
ax2.set_xticklabels(sym_count['Symmetry'], rotation=rot) 
ax2.set(ylim=(0, 100)) 

ax2.grid(axis = 'y', alpha=0.5) 
ax2.set_axisbelow(True) 


ax2.set_ylabel('Frequency Percentage')
#ax2.set_xlabel('Symmetry')

plt.savefig('../v5/sym_frequency_v5.pdf', transparent=False, bbox_inches='tight') 
plt.savefig('../v5/sym_frequency_v5.png', transparent=False, bbox_inches='tight')