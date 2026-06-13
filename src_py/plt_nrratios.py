# Plots g(r) and n(r)
# Import modules
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import warnings
import aux_plt as aux
#------------------------------------------------------------------

# Color/line data; figure defaults
orange = '#FFA500'; dark_g = '#006400'; brown = '#8B4513'
clr_arr = [dark_g,orange,brown,'b','k','m']
mrk_arr = ['s','o','^','d','s','v']
lne_arr = ['-','--']
plt.rc('legend',fontsize=16) # fontsize of the legends
plt.rcParams.update({'font.size': 14}) # other fonts
plt.rcParams.update({'figure.autolayout': True})
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Arial'] + plt.rcParams['font.serif']

# Inputs
dirkey   = 'rdf'
fecharge = 3
trialnum = 5
molality = ['0.1','0.2','0.4','1.0','2.5','5.0']

# RDF inputs
rdf_refatom = 'Fe'
rdf_keys    = {'Fe':1,'NTFSI':2,'OTFSI':3,'FTFSI':4,'OW':7}
rat_cols    = [['OTFSI','OW'],['OTFSI','FTFSI']]
rcol  = 0
ref_rval     = 0.53 # value for ratios are computed

# Directory info
dirsuffix  = 'results_all_trial'+str(trialnum)+'/fetfsi_'+str(fecharge)
headdir    = '../../FeTFSI/' + dirsuffix
resultdir  = headdir + '/' + dirkey + '_all'
figdir     = '../../FeTFSI/figures/' + dirsuffix +'/' + dirkey + '_all'
outheaddir = '../../FeTFSI/analyzed_results/' +  dirsuffix + '/' + dirkey + '_results'

if not os.path.isdir(headdir):
    raise RuntimeError(f'No simulation head directory: {headdir} found')

if not os.path.isdir(figdir):
    os.mkdir(figdir)

if not os.path.isdir(outheaddir):
    os.mkdir(outheaddir)


# Define axes labels for distribution
fig1, ax1 = plt.subplots()
gxlabel = f'Concentration(m)'
gylabel = r'$\phi_{\mathrm{Fe-i}}$($r = r_\mathrm{c,FT}$)'
ax1.set_xlabel(gxlabel,fontsize=16)
ax1.set_ylabel(gylabel,fontsize=16)

for row_index,row in enumerate(rat_cols):

    sel_atom_1 = row[0]
    sel_atom_2 = row[1]

    col_atom_1 = rdf_keys[sel_atom_1]
    col_atom_2 = rdf_keys[sel_atom_2]

    nrat_vals  = np.zeros(len(molality))


    
    for inum,molval in enumerate(molality): # loop in molval

        print(f"Analyzing {molval} m")

        anadir = resultdir + '/rdf_mol_' + str(molval)
        if not os.path.isdir(anadir):
            warnings.warn(f'{anadir} does not exist')
            continue


        # Check file
        nfname  = anadir + '/nrdf_' + rdf_refatom + '_all.xvg'
        if not os.path.exists(nfname):
            warnings.warn(f"{nfname} does not exist in {anadir} ")
            continue

        nralldata = aux.read_xvg_filedata(nfname)
        rindval  = np.abs((nralldata[:,rcol] - ref_rval)).argmin()
        nrat_vals[inum] = nralldata[rindval,col_atom_1]/nralldata[rindval,col_atom_2]


    ax1.plot(np.array(molality,dtype=float),nrat_vals,marker=mrk_arr[row_index],markersize=10,\
             color=clr_arr[row_index],linewidth=2.5,linestyle='--',\
             label=rf'$n_{{\mathrm{{{sel_atom_1}}}}}/n_{{\mathrm{{{sel_atom_2}}}}}$')
ax1.set_ylim(0,2)
ax1.legend()
plt.style.use('seaborn-colorblind')
plt.tight_layout()
fig1.savefig(figdir + '/nofrrat_' + rdf_refatom + '.png',dpi=fig1.dpi)
