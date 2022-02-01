import math
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import scipy.constants as sc
import pandas as pd

path = '/ihome/ghutchison/blp62/GA/running_GA/full_analysis_data.csv'

def split_filename(df):
    df2 = df.copy()

    new = df2['filename'].str.split('_', expand=True)
    df2['acc_term'] = new[0]
    # fixes issue where unit 45 and 90 are the same unit
    #for x in range(len(df2['acc_term'])):
       # if df2['acc_term'][x] == '45':
           # df2['acc_term'][x] = '90'
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
    count_at.to_csv('acc_term_freq_v1.csv')

    count_dc = pd.DataFrame(count_dc.reset_index())
    count_dc.columns = ['don_core', 'don_core_count']
    count_dc.to_csv('don_core_freq_v1.csv')

    count_ac = pd.DataFrame(count_ac.reset_index())
    count_ac.columns = ['acc_core', 'acc_core_count']
    count_ac.to_csv('acc_core_freq_v1.csv')

    return [count_at, count_dc, count_ac]  

run_data = pd.read_csv('/ihome/ghutchison/blp62/GA/running_GA/full_analysis_data.csv', index_col=False, encoding='unicode_escape')

at_counts, dc_counts, ac_counts = get_freq_counts(run_data)

