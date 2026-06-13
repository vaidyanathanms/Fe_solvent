# Import modules
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
#------------------------------------------------------------------

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
fecharge = 2
trialnum = 4
rdfdirs  = [1,2,3]
molality = ['0.1','0.2','0.4','1.0']

# Directories
headdir  = '../../FeTFSI/results_all_trial' + str(trialnum) + '/rdf_fetfsi' + str(fecharge)
figdir   = '../../FeTFSI/figures/rdf_all'
outdir   = '../../FeTFSI/analyzed_results/rdf_all'

refatom  = 'OTFSI'
rcol  = 0
Fecol = 1; rcut_OA_Fe = 0.48
OAcol = 3; rcut_OA_OA = 0.36
HWcol = 8; rcut_OA_HW = 0.26

# Define axes labels for distribution
fig1,ax1 = plt.subplots()
ax1.set_xlabel(r'$r$ (nm)')
ax1.set_ylabel(r'$g_{\mathrm{O(TFSI)-Fe}}(r)$')
plt.style.use('seaborn-colorblind')
plt.tight_layout()

fig2,ax2 = plt.subplots()
ax2.set_xlabel(r'$r$ (nm)')
ax2.set_ylabel(r'$g_{\mathrm{O(TFSI)-O(TFSI)}}(r)$')
plt.style.use('seaborn-colorblind')
plt.tight_layout()

fig3,ax3 = plt.subplots()
ax3.set_xlabel(r'$r$ (nm)')
ax3.set_ylabel(r'$g_{\mathrm{O(TFSI)-H(Water)}}(r)$')
plt.style.use('seaborn-colorblind')
plt.tight_layout()

# coordination arrays

navgOA_Fe = np.zeros(len(molality),dtype=float)
navgOA_OA = np.zeros(len(molality),dtype=float)
navgOA_HW = np.zeros(len(molality),dtype=float)


for inum, molval in enumerate(molality): # loop in runarr

    print(f"Analyzing {molval}")
    nOA_Fe = 0; nOA_OA = 0; nOA_HW = 0
    cntOA_Fe = 0; cntOA_OA = 0; cntOA_HW = 0
    
    for casenum,rdfcase in enumerate(rdfdirs):
    
        anadir = headdir + '/rdfall_mol_' + molval + '/rdfall' + \
            str(rdfcase) + '_' + molval
    
        # Check file
        fname  = anadir + '/nrdf_' + refatom + '_all.xvg'
        if not os.path.exists(fname):
            print(fname, " does not exist! ")
            continue

        # Open and parse file
        with open(fname) as fin:
            lines = (line.lstrip() for line in fin \
                     if not line.lstrip().startswith('#') and \
                     not line.lstrip().startswith('@'))
            data  = np.loadtxt(lines)

            mask_OA_Fe = data[:,rcol] >= rcut_OA_Fe
            mask_OA_OA = data[:,rcol] >= rcut_OA_OA
            mask_OA_HW = data[:,rcol] >= rcut_OA_HW
            
            if np.any(mask_OA_Fe):
                nOA_Fe += data[np.argmax(mask_OA_Fe),Fecol]
                cntOA_Fe += 1
            if np.any(mask_OA_OA):
                nOA_OA += data[np.argmax(mask_OA_OA),OAcol]
                cntOA_OA += 1
            if np.any(mask_OA_HW):
                nOA_HW += data[np.argmax(mask_OA_HW),HWcol]
                cntOA_HW += 1
                        
    # Divide ydata by rdata
    navgOA_Fe[inum] = nOA_Fe/cntOA_Fe
    navgOA_OA[inum] = nOA_OA/cntOA_OA
    navgOA_HW[inum] = nOA_HW/cntOA_HW

# Bar and figure settings    
width = 0.25
x = np.arange(len(molality))
fig, ax = plt.subplots(figsize=(8, 5))


# Plot bars
plt.bar(x - width, navgOA_Fe, width, label=r'O$_{TFSI}$-Fe' + '; r$_{cut}$: ' + str(rcut_OA_Fe) + ' nm')
plt.bar(x,         navgOA_OA, width, label=r'O$_{TFSI}$-O$_{TFSI}$' +  '; r$_{cut}$: ' + str(rcut_OA_OA) + ' nm')
plt.bar(x + width, navgOA_HW, width, label=r'O$_{TFSI}$-H$_{Water}$'+  '; r$_{cut}$: ' + str(rcut_OA_Fe) + ' nm')

# Formatting
plt.xticks(x, molality)
plt.xlabel("Molality")
plt.ylabel("Coordination Number")
# Legend INSIDE top-right
leg = ax.legend(
    loc='upper right',
    frameon=True,
    framealpha=0.9,
    borderaxespad=0.4,   # small gap to axes edges
    labelspacing=0.4     # tighter label spacing for long legends
)
# --- Auto add extra vertical margin to avoid overlap with tall legend ---
fig.canvas.draw()  # needed to get accurate legend size
bbox_display = leg.get_window_extent(fig.canvas.get_renderer())

# Convert legend bbox to axes (0..1) coordinates
bbox_axes = bbox_display.transformed(ax.transAxes.inverted())
legend_height_frac = bbox_axes.height  # fraction of axes height occupied by legend

# Add padding proportional to legend height
pad = legend_height_frac + 0.02  # a little extra
ymax = np.max(np.maximum.reduce([navgOA_Fe,navgOA_OA,navgOA_HW]))
ax.set_ylim(0, 1.1*ymax * (1 + pad))

plt.tight_layout()
plt.savefig('../../FeTFSI/figures/rdf_all/OA_coordnum_Fe' + str(fecharge) + '.png')

# Write to output
output_data = np.column_stack((np.asarray(molality,dtype=str), navgOA_Fe, navgOA_OA, navgOA_HW)) # converts to str!
np.savetxt(outdir + '/CN_' + refatom + 'Fe_' + str(fecharge) +  '.txt', output_data,\
           header="Molality   navgOA_Fe   navgOA_OA   navgOA_HW",\
           fmt=['%s', '%s', '%s', '%s'], delimiter = '\t',comments='')

