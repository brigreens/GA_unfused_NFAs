import csv
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

SMALL_SIZE = 10
MEDIUM_SIZE = 14 
BIGGER_SIZE = 18 
plt.rc('font', size=BIGGER_SIZE)          # controls default text sizes 
plt.rc('axes', titlesize=BIGGER_SIZE)     # fontsize of the axes title 
plt.rc('axes', labelsize=BIGGER_SIZE)    # fontsize of the x and y labels 
plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels 
plt.rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels 
plt.rc('legend', fontsize=MEDIUM_SIZE)    # legend fontsize 
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title"


df_v1 = pd.read_csv('all_NFAs_PCE_v1.csv')
df_v2 = pd.read_csv('all_NFAs_PCE_v2.csv')
df_v3 = pd.read_csv('all_NFAs_PCE_v3.csv')
df_v4 = pd.read_csv('all_NFAs_PCE_v4.csv')
df_all = pd.read_csv('all_NFAs_PCE_v1245.csv')


'''
plt.hist(df_v1['PCE'], bins = 100, color = 'r', alpha = 0.25, label = 'Version 1')
plt.hist(df_v2['PCE'], bins = 100, color = 'orange', alpha = 0.25, label = 'Version 2')
plt.hist(df_v3['PCE'], bins = 100, color = 'green', alpha = 0.25, label = 'Version 3')
plt.hist(df_v4['PCE'], bins = 100 , color = 'blue', alpha = 0.25, label = 'Version 4')
#plt.hist(df_all['PCE'], bins = 100, color = 'b', alpha = 0.5)
plt.legend(loc = 'upper left')
plt.xlabel('PCE')
plt.ylabel('Frequency of Occurence')
plt.savefig('hist_NFAs_all_oneplot.png')
plt.clf()'''


sns.distplot(df_v1, hist=False, kde=True, kde_kws = {'linewidth': 3}, label = 'Version 1')
sns.distplot(df_v2, hist=False, kde=True, kde_kws = {'linewidth': 3}, label = 'Version 2')
sns.distplot(df_v3, hist=False, kde=True, kde_kws = {'linewidth': 3}, label = 'Version 3')
sns.distplot(df_v4, hist=False, kde=True, kde_kws = {'linewidth': 3}, label = 'Version 4')
plt.axvline(x=18, linestyle = '--', color = 'black', alpha = 0.5)

plt.legend()
plt.xlabel('PCE (%)')
plt.ylabel('Density')
plt.tight_layout
plt.savefig('KDE_allversions.png', bbox_inches = "tight")
plt.savefig('KDE_allversions.pdf', dpi=600, bbox_inches = "tight")

plt.clf()

fig, ax = plt.subplots()
N, bins, patches = ax.hist(df_all['PCE'], bins = 100)

for i in range(0, 74):
    patches[i].set_facecolor('tab:red')
    patches[i].set_alpha(0.7)

for i in range(74, len(patches)):
    patches[i].set_facecolor('tab:green')
    patches[i].set_alpha(0.7)

plt.xlabel('PCE')
plt.ylabel('Frequency of Occurence')
plt.savefig('hist_NFAs_all.png', bbox_inches = "tight")
plt.savefig('hist_NFAs_all.pdf', dpi=600, bbox_inches = "tight")
plt.clf()

'''plt.hist(df_all['below18'], bins = 94, color = 'r')
plt.hist(df_all['above18'], bins = 26, color = 'g')
plt.xlabel('PCE')
plt.ylabel('Frequency of Occurence')
plt.savefig('hist_NFAs_all.png')
plt.savefig('hist_NFAs_all.pdf')
plt.clf()'''

'''
plt.xlabel('PCE')
plt.ylabel('Frequency of Occurence')
plt.savefig('hist_NFAs_v2.png')
plt.savefig('hist_NFAs_v2.pdf')
plt.clf()

plt.xlabel('PCE')
plt.ylabel('Frequency of Occurence')
plt.savefig('hist_NFAs_v3.png')
plt.savefig('hist_NFAs_v3.pdf')
plt.clf()

plt.xlabel('PCE')
plt.ylabel('Frequency of Occurence')
plt.savefig('hist_NFAs_v4.png')
plt.savefig('hist_NFAs_v4.pdf')
plt.clf()

plt.savefig('hist_NFAs_all.png')
plt.savefig('hist_NFAs_all.pdf')

fig, ax = plt.subplots(dpi=300) '''
