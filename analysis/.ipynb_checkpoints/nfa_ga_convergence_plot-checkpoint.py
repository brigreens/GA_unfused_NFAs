import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams


df_max = pd.read_csv("../opt_GA/quick_analysis_data.csv", index_col=False)

SMALL_SIZE = 8
MEDIUM_SIZE = 10
BIGGER_SIZE = 12

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=MEDIUM_SIZE)  # fontsize of the figure title
rcParams['text.usetex'] = True

# PCE
min_pce_m = min(df_max[' max_PCE'])
max_pce_m = max(df_max[' max_PCE'])

# deltaHOMO
min_homo_m = min(df_max['max_deltaHOMO'])
max_homo_m = max(df_max['max_deltaHOMO'])

# summedoscs
min_sum_m = min(df_max['max_summedoscs'])
max_sum_m = max(df_max['max_summedoscs'])

print('min_pce_term', min_pce_m)
print('max_pce_term', max_pce_m)
print('min_homo_term', min_homo_m)
print('max_homo_term', max_homo_m)
print('min_sum_term', min_sum_m)
print('max_sum_term', max_sum_m)

fig, ((ax1), (ax2), (ax3)) = plt.subplots(3, 1, sharex=True, dpi=300)
#fig, (ax1, ax2, ax3) = plt.subplots(1, 3, dpi=300)


# set figure axes
xlimits = (0, 115)
xticks = []
xticks_gen = [0, 30, 60, 90, 115]

# PCE y axis
ylim_p = (15, 25)
yticks_p = [15, 20, 25]

# deltaHOMO y axis
ylim_h = (0, 2)
yticks_h =[0, 0.5, 1, 1.5]

# summedoscs y axis
ylim_s = (6, 12)
yticks_s =[6, 8, 10, 12]

ax1.plot('gen', ' max_PCE', data=df_max, linestyle='', marker='o', markersize=2)
ax1.set_ylabel('Best PCE', fontsize=8)
ax1.set(ylim=ylim_p, xlim=xlimits)
ax1.set_yticks(yticks_p)
ax1.set_xlabel('Generation', fontsize=8)


ax1.set_xticks(xticks)
ax1.grid(alpha=0.5)
ax1.set_axisbelow(True)

#for tic in ax1.xaxis.get_major_ticks():
    #tic.tick1On = tic.tick2On = False

ax2.plot('gen', 'max_deltaHOMO', data=df_max, linestyle='', marker='o', markersize=2,color='#cb2121')
ax2.set_ylabel('Best $\Delta$HOMO', fontsize=8)
ax2.set(ylim=ylim_h, xlim=xlimits)
ax2.set_yticks(yticks_h)
ax2.set_xlabel('Generation', fontsize=8)


ax2.set_xticks(xticks)
ax2.grid(alpha=0.5)
ax2.set_axisbelow(True)

#for tic in ax2.xaxis.get_major_ticks():
    #tic.tick1On = tic.tick2On = False

ax3.plot('gen', 'max_summedoscs', data=df_max, linestyle='', marker='o', markersize=2,color='#cb2121')
ax3.set_ylabel('Best $\sum{f}$', fontsize=8)
ax3.set(ylim=ylim_s, xlim=xlimits)
ax3.set_yticks(yticks_s)
ax3.set_xlabel('Generation', fontsize=8)


ax3.set_xticks(xticks_gen)
ax3.grid(alpha=0.5)
ax3.set_axisbelow(True)

    
top = 0.75
bottom = 0.12
midpoint = (top+bottom)/2

left = 0.11
right = 0.35
    
#fig.text((left+right)/2, 0.02,'Generation', ha='center')

plt.subplots_adjust(top=top, bottom=bottom, left=left, right=right, hspace=0.3)

plt.savefig('../opt_GA/convergence_plots.pdf', transparent=False, bbox_inches='tight')
plt.savefig('../opt_GA/convergence_plots.png', transparent=False, bbox_inches='tight')

