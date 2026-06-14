# Plots S(q) and partial S(q)
# How to use
# 1. Use the scripts in src_mat to compute S(q) and partial S(q)
# 2. Check whether the files are saved in
# ../../FeTFSI/analyzed_results/sofq_all


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
mrk_arr = ['x','o','^','d','s','v']
lne_arr = ['-','--']
plt.rc('legend',fontsize=16) # fontsize of the legends
plt.rcParams.update({'font.size': 14}) # other fonts
plt.rcParams.update({'figure.autolayout': True})
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Arial'] + plt.rcParams['font.serif']

# Inputs
dirkey   = 'sofq'
fecharge = 3
trialnum = 5
molality = ['0.1','0.2','0.4','1.0','2.5','5.0']

# Sofq inputs
sofq_pairs = ['FF','AA','AF','WF','WA','WW']
maxylimplt = [5,5,15,15,10]

# Directory info
dirsuffix  = 'analyzed_results/results_all_trial'+str(trialnum)+'/fetfsi_'+str(fecharge)
headdir    = '../../FeTFSI/' +  dirsuffix
resultdir  = headdir + '/' + dirkey + '_all'
figdir     = '../../FeTFSI/figures/' + dirsuffix +'/' + dirkey + '_all'

# Define axes labels for distribution
for inum,molval in enumerate(molality): # loop in molval

    print(f"Analyzing {molval} m")    
    # Plot S(q) data
    fig1, ax1 = plt.subplots()
    gxlabel = r'$q$ (nm$^{-1}$)'
    gylabel = r'$I(q)$'
    ax1.set_xlabel(gxlabel,fontsize=16)
    ax1.set_ylabel(gylabel,fontsize=16)
    plt.style.use('seaborn-colorblind')
    plt.tight_layout()

    # Check file
    Itotfname   = anadir + '/avgsofq_' + rdf_refatom + '_all.xvg'
    if not os.path.exists(gfname):
        warnings.warn(f"{gfname} does not exist in {anadir} ")
        continue

         
    # Plot partial S(q) data
    fig2, ax2 = plt.subplots()
    gxlabel = r'$q$ (nm$^{-1}$)'
    gylabel = r'$I_{AB}(q)$'
    ax2.set_xlabel(nxlabel,fontsize=20)
    ax2.set_ylabel(nylabel,fontsize=20)
    plt.style.use('seaborn-colorblind')
    plt.tight_layout()

    for sid,spair in enumerate(sofq_pairs):

        print(f'Analyzing {spair}')

        # Check file
        Ipartfname  = anadir + '/psofq_' + rdf_refatom + '_all.xvg'
        if not os.path.exists(nfname):
            warnings.warn(f"{nfname} does not exist in {anadir} ")
            continue

        # Open and parse g(r) file
        gralldata = aux.read_xvg_filedata(gfname)
        # Open and parse n(r) file
        nralldata = aux.read_xvg_filedata(nfname)
            
                
        rlendata = len(gralldata[:,rcol])
        rxgrdata = np.zeros(rlendata)
        gofrdata = np.zeros(rlendata)

        rlendata = len(nralldata[:,rcol])
        rxnrdata = np.zeros(rlendata)
        nofrdata = np.zeros(rlendata)

        rxgrdata += gralldata[:,rcol]
        gofrdata += gralldata[:,rdfcolval]
        rxnrdata += nralldata[:,rcol]
        nofrdata += nralldata[:,rdfcolval]
        #nfylcnt += 1

        # Divide ydata by rdata
        #rxgrdata   /= nfylcnt; gofrdata /= nfylcnt; rxnrdata /= nfylcnt;  nofrdata /= nfylcnt

        if rdfkey == 'Fe':
            ax1.plot(rxgrdata,gofrdata,color=clr_arr[inum],marker=mrk_arr[inum],linestyle='None',label=molval + ' m')
        else:
            ax1.plot(rxgrdata,gofrdata,color=clr_arr[inum],linestyle='-',label=molval + ' m')
        ax2.plot(rxnrdata,nofrdata,color=clr_arr[inum],linestyle='-', linewidth=2.5,label=molval + ' m')
        ax3.plot(rxgrdata,gofrdata,color=clr_arr[inum],linestyle='-', linewidth=2.5,label=molval + ' m')
        ax3.plot(rxnrdata,nofrdata,color=clr_arr[inum],linestyle='--',linewidth=2.5,label='_Hidden')

    ax1.legend()
    ax2.legend()
    ax3.legend()
    
    ax1.set_xlim(0,maxxlimplt[cid])
    ax2.set_xlim(0,maxxlimplt[cid]); ax2.set_ylim(0,maxylimplt[cid])
    ax3.set_xlim(0,maxxlimplt[cid]); ax3.set_ylim(0,maxylimplt[cid])

    fig1.savefig(figdir + '/gofr_' + rdf_refatom + '_' + rdfkey + '.png',dpi=fig1.dpi)
    fig2.savefig(figdir + '/nofr_' + rdf_refatom + '_' + rdfkey + '.png',dpi=fig2.dpi)
    fig3.savefig(figdir + '/combgnofr_' + rdf_refatom + '_' + rdfkey + '.png',dpi=fig3.dpi)
    
    plt.close(fig1)
    plt.close(fig2)
    plt.close(fig3)



