import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

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

# Bar width and x positions
bar_width = 0.35

df = pd.read_excel('../../FeTFSI/Diffusivity/diff_calc.xlsx', sheet_name='Diff_Fatoms')

mol_data = df['molality']
D_Fe2    = df['D_Fe(TFSI)2']
Ebar_Fe2 = df['Ebar_Fe(TFSI)2']
D_Fe3    = df['D_Fe(TFSI)3']
Ebar_Fe3 = df['Ebar_Fe(TFSI)3']

x = np.arange(len(mol_data))

# Create the grouped bar plot
h = plt.figure(figsize=(8, 6))
plt.bar(x - bar_width/2, D_Fe2, bar_width, yerr=Ebar_Fe2, capsize=5, label=r'Fe(TFSI)$_2$')
plt.bar(x + bar_width/2, D_Fe3, bar_width, yerr=Ebar_Fe3, capsize=5, label=r'Fe(TFSI)$_3$')

# Formatting
plt.xlabel('Molality')
plt.ylabel(r'Diffusivity $\times 10^{10}$ (m$^2$/s)')
plt.xticks(ticks=x, labels=mol_data)
plt.legend()
plt.tight_layout()
plt.savefig('../../FeTFSI/figures/diff_all/diff_fatoms.png')
