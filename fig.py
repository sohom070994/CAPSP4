# -*- coding: utf-8 -*-
"""
Created on Mon May  1 10:53:31 2023

@author: hebaish
"""

#%%
import matplotlib.pyplot as plt
import os
import numpy as np
import pandas as pd
import scipy.stats as st
from scipy import stats



dir_path = os.path.dirname(os.path.realpath(__file__))
os.chdir(dir_path)
plt.rcParams['font.family'] = "Times New Roman"

# Generate some data
#%%
# Import data
file_names = ['0mod_Model_E3_0mod_ResponseDetail.csv',
              '1mod_Model_E3_1mod_ResponseDetail.csv', 
              '2mod_Model_E3_2mod_ResponseDetail.csv', 
              '3mod_Model_E3_3mod_ResponseDetail.csv', 
              '4mod_Model_E3_4mod_ResponseDetail.csv', 
              '5mod_Model_E3_5mod_ResponseDetail.csv', 
              '6mod_Model_E3_6mod_ResponseDetail.csv']

dataframes = {}
# Loop through the file names
opt_triage = 33.40
opt_crisis = 0.815
uni_triage = 71.10
uni_crisis = 1.15
triage_means = []
crisis_means = []
triage_CI = []
crisis_CI = []
wtd_means = []
wtd_CI = []
improv_means = []
improv_CI = []


mod=0
mods = []
for file_name in file_names:
    mods.append(mod)
    # Read the CSV file into a DataFrame
    df = pd.read_csv(file_name)
    
    triage_means.append(df['RAccessTime'].mean()/8)
    crisis_means.append(df['CAccess'].mean())
    triage_CI.append((df['RAccessTime'].mean()-st.t.interval(alpha=0.95, df=len(df['RAccessTime'])-1, loc=np.mean(df['RAccessTime']), scale=st.sem(df['RAccessTime']))[0])/8)
    crisis_CI.append(df['CAccess'].mean()-st.t.interval(alpha=0.95, df=len(df['CAccess'])-1, loc=np.mean(df['CAccess']), scale=st.sem(df['CAccess']))[0])
    
    if mod > 0:
        wtd_means.append(df['RAccessTime'].mean()/opt_triage + df['CAccess'].mean()/opt_crisis)
        wtd = df['RAccessTime']/opt_triage + df['CAccess']/opt_crisis
        wtd_CI.append(np.mean(wtd)-st.t.interval(alpha=0.95, df=len(wtd)-1, loc=np.mean(wtd), scale=st.sem(wtd))[0])
        
        triage_improv = abs(df['RAccessTime'].mean()-uni_triage)/uni_triage
        crisis_improv = abs(df['CAccess'].mean()-uni_crisis)/uni_crisis
        improv_means.append((triage_improv+crisis_improv)/2)
        improv = (abs(df['RAccessTime']-uni_triage)/uni_triage + abs(df['CAccess']-uni_crisis)/uni_crisis)/2
        improv_CI.append(np.mean(improv)-st.t.interval(alpha=0.95, df=len(improv)-1, loc=np.mean(improv), scale=st.sem(improv))[0])
        
    # Store the DataFrame in a separate variable with the corresponding name
    #globals()['df' + str(df_number)] = df
    dataframes['df' + str(mod)] = df
    mod+=1
#%%
x = [0,1,2,3,4,5,6]

y1 = [8.8875, 4.9875,4.925,4.45,4.575,4.2125,4.2625]

y2 = [1.15,0.793,0.784,0.85,0.82,0.82,0.79]

ci1 = [0.53655,0.2793,0.2793,0.239365,0.25725,0.204085,0.197225]

ci2 = [0.14308,0.0882,0.0882,0.10388,0.09996,0.098,0.08036]


# Create the figure
fig, ax1 = plt.subplots(figsize=(8, 5), dpi=300)
ax2 = ax1.twinx()

# Plot the data with error bars
ax1.plot(x, y1, color='black', label='Triage')
ax2.plot(x, y2, color='grey', ls='--', label='Crisis')

ax1.errorbar(x, y1, yerr=ci1, color='black', capsize=4,marker='o')
ax2.errorbar(x, y2, yerr=ci2, color='grey', capsize=4,marker='o', ls='--')

# Set axis labels and font
ax1.set_xlabel('Allowable modifications ($M$)', fontsize=20, fontfamily="Times New Roman")
ax1.set_ylabel('Expected triage access time (days)', fontsize=20, color='black')
ax2.set_ylabel('Expected crisis access time (hours)', fontsize=20, color='black')
ax1.tick_params(axis='y', labelcolor='black')
ax2.tick_params(axis='y', labelcolor='black')
ax1.set_ylim(3.5, 10)
ax2.set_ylim(0, 1.8)
for ax in (ax1, ax2):
    for tick in ax.get_xticklabels() + ax.get_yticklabels():
        tick.set_fontname('Times New Roman')
        tick.set_fontsize(18)

# Add legend
lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
lines = lines1 + lines2
labels = labels1 + labels2
ax1.legend(lines, labels, loc='best', prop={'family': 'Times New Roman', 'size': 18})

# Save the figure in pdf format
#plt.savefig('fig8b.pdf', format='pdf', dpi=300, bbox_inches='tight')

# %%
x2= [1,2,3,4,5,6]
y3= [2.16761691341244,2.14160390874692,2.10881304874913,2.10194335255869,2.0151170052533,1.99028323720657]
y4= [0.387641726297655,0.394990478068423,0.404253994281652,0.406194704041762,0.430723407346676,0.437739021238837]
# Create the figure
fig2, ax3 = plt.subplots(figsize=(8, 5), dpi=300)
ax4 = ax3.twinx()

# Plot the data with error bars
ax3.plot(x2, y3, color='black', label='Weighted access time', alpha=0.7)
ax4.plot(x2, [x*100 for x in y4], color='grey', ls='--', label='% improvement', alpha=1)

ax3.errorbar(x2, y3, yerr=wtd_CI, color='black', capsize=4,marker='o', alpha=1)
ax4.errorbar(x2, [x*100 for x in y4], yerr=[x*100 for x in improv_CI], color='grey', capsize=4,marker='o', ls='--', alpha=1)


# Set axis labels and font
ax3.set_xlabel('Allowable modifications ($M$)', fontsize=20)
ax3.set_ylabel('Weighted access time', fontsize=20, color='black')
ax4.set_ylabel('Percentage improvment (%)', fontsize=20, color='black')
from matplotlib.ticker import FormatStrFormatter

ax3.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))

ax3.tick_params(axis='y', labelcolor='black')
ax4.tick_params(axis='y', labelcolor='black')
ax3.set_ylim(1.9, 2.25)
ax4.set_ylim(33, 48)
for ax in (ax3, ax4):
    for tick in ax.get_xticklabels() + ax.get_yticklabels():
        tick.set_fontname('Times New Roman')
        tick.set_fontsize(18)

# Add legend
lines3, labels3 = ax3.get_legend_handles_labels()
lines4, labels4 = ax4.get_legend_handles_labels()
lines2 = lines3 + lines4
labels2 = labels3 + labels4
ax3.legend(lines2, labels2, loc='lower center', prop={'family': 'Times New Roman', 'size': 18})
# plt.savefig('fig8a.pdf', format='pdf', dpi=300, bbox_inches='tight')

# %%
