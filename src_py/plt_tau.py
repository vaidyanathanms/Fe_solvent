# To fit tau for autocorrelation functions
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import glob
import warnings
from pathlib import Path
import aux_plt as aux
import math

system = 'cades' # cades or lap

# Color/line data; figure defaults

orange = '#FFA500'; dark_g = '#006400'; brown = '#8B4513'
clr_arr = [dark_g,orange,brown,'b','k','m']
mrk_arr = ['o','d','s','h','x','^','p']
lne_arr = ['-','--']
plt.rc('legend',fontsize=16) # fontsize of the legends
plt.rcParams.update({'font.size': 14}) # other fonts
plt.rcParams.update({'figure.autolayout': True})
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Arial'] + plt.rcParams['font.serif']

# Inputs
fecharge = 3
trialnum = 5
molality = ['0.1','0.2','0.4','1.0','2.5']
iontype  = ['ion']

# Directory details
dirkey     = 'tau'
dirsuffix  = 'results_all_trial'+str(trialnum)+'/fetfsi_'+str(fecharge)
headdir    = '../../FeTFSI/' + dirsuffix
figdir     = '../../FeTFSI/figures/' + dirsuffix +'/' + dirkey + '_all'
anadir = '../../FeTFSI/analyzed_results/' + dirsuffix + '/tau_results'


if not os.path.isdir(figdir):
    os.mkdir(figdir)


if not os.path.isdir(anadir):
    raise RuntimeError(f'{anadir} does not exist')


fig1,ax1 = plt.subplots()
ax1.set_xlabel(r'$q$ ($\AA^{-1}$)',fontsize=16)
ax1.set_ylabel(r'$\tau_{c}$ (ps)',fontsize=16)
plt.style.use('seaborn-colorblind')
plt.tight_layout()

fig2,ax2 = plt.subplots()
ax2.set_xlabel(r'$q$ ($\AA^{-1}$)',fontsize=16)
ax2.set_ylabel(r'$\beta$',fontsize=16)
plt.style.use('seaborn-colorblind')
plt.tight_layout()


# Load data
for atype in iontype:
    
    for inum,molval in enumerate(molality): # loop in runarr
       
        print(f"Analyzing {atype} tau for {molval}")

        if float(molval) > 1.0:
            rtype = 'main'
        else:
            rtype = 'smallrun'
            
        # Check file
        list_fnames = glob.glob(f'{anadir}/tau{atype}*{rtype}*{molval}*')
        if list_fnames == []:
            warnings.warn(f'No tau{atype}*{rtype}* files exist in {anadir}')
            continue
        fname = max(list_fnames,key=os.path.getmtime)

        data = np.loadtxt(fname,comments = "#")

        qval = data[:,0]
        beta = data[:,2]
        tau_avg = data[:,3]
        tau_err = data[:,4]

        mask = np.isfinite(qval) & np.isfinite(tau_avg) & np.isfinite(tau_err)
        
        qval = qval[mask]
        beta = beta[mask]
        tau_avg = tau_avg[mask]
        tau_err = tau_err[mask]

        ax1.errorbar(qval,tau_avg,yerr=tau_err,marker=mrk_arr[inum],\
                    markersize=6,linestyle='None',capsize=4,\
                    color=clr_arr[inum],label=f'{molality[inum]} m')

        ax2.plot(qval,beta,marker=mrk_arr[inum],\
                 markersize=6,linestyle='None',\
                 color=clr_arr[inum],label=f'{molality[inum]} m')


ax1.legend(loc = 'lower left')
ax2.legend()
ax1.set_yscale('log')
# Set output file and output figure data
fig1.savefig(f'{figdir}/tau_{atype}_{rtype}.png',dpi=fig1.dpi)
fig2.savefig(f'{figdir}/beta_{atype}_{rtype}.png',dpi=fig2.dpi)
