import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import warnings


# Color/line data; figure defaults
orange = '#FFA500'; dark_g = '#006400'; brown = '#8B4513'
clr_arr = [dark_g,orange,brown,'b','k','m']
mrk_arr = ['o','d','s']
lne_arr = ['-','--']
plt.rc('legend',fontsize=16) # fontsize of the legends
plt.rcParams.update({'font.size': 14}) # other fonts
plt.rcParams.update({'figure.autolayout': True})
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']

# Inputs
fecharge = 3
trialnum = 5

# Read the Excel file
inpfile = '../../FeTFSI/analyzed_results/density_all/density.xlsx'
df = pd.read_excel(inpfile,sheet_name='Fe'+str(fecharge)+
                   '_trial_'+str(trialnum))

# print column names
print(df.columns)

# Extract columns
x = df["mol"]
y_fe = df["d_Fe(cm2/s)"]
err_fe = df["ed_Fe(cm2/s)"]
y_f = df["d_F(cm2/s)"]
err_f = df["ed_F(cm2/s)"]

figdir   = '../../FeTFSI/figures/diff_all'

if not os.path.isdir(figdir):
    os.mkdir(figdir)


# Create scatter plot with error bars
fig1,ax1 = plt.subplots()
ax1.set_xlabel(r'Concentration (m)')
ax1.set_ylabel(r'Diffusivity (m$^2$/s)')
plt.style.use('seaborn-colorblind')
plt.tight_layout()

plt.errorbar(
    x, 0.0001*y_fe, yerr=0.0001*err_fe,
    fmt='o', capsize=4, linestyle='none',
    label='Fe'
)

plt.errorbar(
    x, 0.0001*y_f, yerr=0.0001*err_f,
    fmt='s', capsize=4, linestyle='none',
    label='F (TFSI)'
)

plt.legend()
fig1.savefig(figdir + '/diff_all_Fe' + str(fecharge) +
             '_trial_' + str(trialnum) + '.png',dpi=fig1.dpi)



# Create scatter plot with error bars in semilog scale
fig2,ax2 = plt.subplots()
ax2.set_xlabel(r'Concentration (m)')
ax2.set_ylabel(r'Diffusivity (m$^2$/s)')
plt.style.use('seaborn-colorblind')
plt.tight_layout()

plt.errorbar(
    x, 0.0001*y_fe, yerr=0.0001*err_fe,
    fmt='o', capsize=4, linestyle='none',
    label='Fe'
)

plt.errorbar(
    x, 0.0001*y_f, yerr=0.0001*err_f,
    fmt='s', capsize=4, linestyle='none',
    label='F (TFSI)'
)

plt.yscale('log')
plt.legend()
fig2.savefig(figdir + '/logdiff_all_Fe' + str(fecharge) +
             '_trial_' + str(trialnum) + '.png',dpi=fig1.dpi)

