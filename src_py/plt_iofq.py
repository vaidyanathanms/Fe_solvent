# Plots I(q) and partial I(q)
# How to use
# 1. Use the scripts in src_mat (compute_Iofq.m) to compute I(q) and partial I(q)
# 2. Check whether the files are saved in
# ../../FeTFSI/analyzed_results/Iofq_all


# Import modules
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd
import os
import warnings
import aux_plt as aux
#------------------------------------------------------------------

# Color/line data; figure defaults
orange = '#FFA500'; dark_g = '#006400'; brown = '#8B4513'
clr_arr = [dark_g,orange,brown,'b','k','m']
mrk_arr = ['x','o','^','d','s','v']
lne_arr = ['-','--']
plt.rc('legend',fontsize=16) # fontsize of the legends
plt.rcParams.update({'font.size': 14}) # other fonts
plt.rcParams.update({'figure.autolayout': True})
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Arial'] + plt.rcParams['font.serif']

# Inputs
dirkey   = 'Iofq'
fecharge = 3
trialnum = 5
molality = ['0.1','0.2','0.4','1.0','2.5','5.0']
dirval = 1
lorwind = 1
normtyp = 'add'

# Iofq inputs
iofq_pairs  = ['FF','AF','WF','AA','WA','WW']
iofq_labels = ['Fe-Fe','Fe-TFSI','Fe-Water','TFSI-TFSI','TFSI-Water','Water-Water']

# Directory info
dirsuffix  = 'results_all_trial'+str(trialnum)+'/fetfsi_'+str(fecharge)
anadir     = '../../FeTFSI/analyzed_results/' +  dirsuffix + '/Iofq_results'
figdir     = '../../FeTFSI/figures/' + dirsuffix +'/' + dirkey + '_all'

if not os.path.isdir(anadir):
    raise RuntimeError(f'{anadir} not found')

# Define axes labels for distribution
# Plot S(q)-all data
fig1, ax1 = plt.subplots()
xlabel = r'$q$ ($\mathrm{\AA}^{-1}$)'
xlabel = r'$q$ (nm$^{-1}$)'
ylabel = r'$I(q)$'
ax1.set_xlabel(xlabel,fontsize=16)
ax1.set_ylabel(ylabel,fontsize=16)
plt.style.use('seaborn-colorblind')
plt.tight_layout()

# Plot S(q)-high data
fig3, ax3 = plt.subplots()
xlabel = r'$q$ ($\mathrm{\AA}^{-1}$)'
xlabel = r'$q$ (nm$^{-1}$)'
ylabel = r'$I(q)$'
ax3.set_xlabel(xlabel,fontsize=16)
ax3.set_ylabel(ylabel,fontsize=16)
plt.style.use('seaborn-colorblind')
plt.tight_layout()


for inum,molval in enumerate(molality): # loop in molval

    print(f"Analyzing {molval} m")    
    # Check file
    Itotfname   = anadir + '/avgIofq_fe_' + str(fecharge) + '_mol_' + str(molval) + \
        '_lorwind_' + str(lorwind) + '_norm_' +  normtyp + '.dat'

    if not os.path.exists(Itotfname):
        warnings.warn(f"{Itotfname} does not exist in {anadir} ")
        continue

    df_itot =  pd.read_csv(Itotfname,sep=r'\s+')
    ax1.plot(df_itot["q"].to_numpy(),df_itot["AvgIofq"].to_numpy(),color=clr_arr[inum],\
             linestyle='-',linewidth=2.5,label=molval + ' m')

    if float(molval) == 1.0 or float(molval) == 2.5:
        ax3.plot(df_itot["q"].to_numpy(),df_itot["AvgIofq"].to_numpy(),color=clr_arr[inum],\
                 linestyle='-',linewidth=2.5,label=molval + ' m')
        
    
    # Plot partial S(q) data
    fig2, ax2 = plt.subplots()
    xlabel = r'$q$ ($\mathrm{\AA}^{-1}$)'
    ylabel = r'$I_{AB}(q)$'
    ax2.set_xlabel(xlabel,fontsize=16)
    ax2.set_ylabel(ylabel,fontsize=16)
    plt.style.use('seaborn-colorblind')
    plt.tight_layout()

    if inum == 0: ymin = df_itot["AvgIofq"].min()
    if inum == len(molality)-1: ymax = df_itot["AvgIofq"].max()
    
    for sid,spair in enumerate(Iofq_pairs):

        print(f'Analyzing {spair}')

        # Check file
        Ipartfname  = anadir + '/pIofq_' + spair + '_Fe_' + str(fecharge) + '_mol_' + str(molval)\
            + '_dir_' + str(dirval) + '_lorwind_' + str(lorwind) + '_norm_' +  normtyp + '.dat'
        if not os.path.exists(Ipartfname):
            warnings.warn(f"{Ipartfname} does not exist in {anadir} ")
            continue

        np_ptot =  np.loadtxt(Ipartfname,skiprows=1,dtype=float)
        ax2.plot(np_ptot[:,0],np_ptot[:,1],color=clr_arr[sid],\
                 linestyle='--', linewidth=2.5,label=Iofq_labels[sid])
                


    ax2.legend()
    ax2.set_xlim(0.25,1.2)
    ax2.set_xscale('log')
    ax2.set_xticks([0.3, 0.4, 0.6, 0.8, 1.0])
    ax2.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: f'{x:.1f}'))
    fig2.savefig(f'{figdir}/ipart_{molval}.png',dpi=fig2.dpi)
    plt.close(fig2)
    
ax1.legend(loc='upper left')
ax1.set_xlim(0.25,1.3)
ax1.set_ylim(ymin-2,ymax+40)
ax1.set_xscale('log')
ax1.set_xticks([0.3, 0.4, 0.6, 0.8, 1.0])
ax1.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: f'{x:.1f}'))
fig1.savefig(f'{figdir}/itot_all.png',dpi=fig1.dpi)
plt.close(fig1)

ax3.legend()
ax3.set_xlim(0.3,1.0)
ax3.set_xscale('log')
ax3.set_xticks([0.3, 0.4, 0.6, 0.8, 1.0])
ax3.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: f'{x:.1f}'))
fig3.savefig(f'{figdir}/itot_highconc.png',dpi=fig1.dpi)
plt.close(fig3)




