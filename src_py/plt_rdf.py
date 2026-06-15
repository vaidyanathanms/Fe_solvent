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
mrk_arr = ['x','o','^','d','s','v']
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
maxxlimplt  = [3.0,2,1.2,1.0,0.3]
maxylimplt = [5,5,15,15,10]
rcol  = 0

# Directory info
dirsuffix  = 'results_all_trial'+str(trialnum)+'/fetfsi_'+str(fecharge)
headdir    = '../../FeTFSI/' + dirsuffix
resultdir  = headdir + '/' + dirkey + '_all'
figdir     = '../../FeTFSI/figures/' + dirsuffix +'/' + dirkey + '_all'

if not os.path.isdir(headdir):
    raise RuntimeError(f'No simulation head directory: {headdir} found')

if not os.path.isdir(figdir):
    os.mkdir(figdir)


# Define axes labels for distribution
for cid, (rdfkey,rdfcolval) in enumerate(rdf_keys.items()):

    print(f'Analyzing {rdf_refatom} - {rdfkey}')
    
    # Plot g(r) data
    fig1, ax1 = plt.subplots()
    gxlabel = rf'$r$ (nm)'
    gylabel = rf'$g_{{\mathrm{{{rdf_refatom}-{rdfkey}}}}}(r)$'
    ax1.set_xlabel(gxlabel,fontsize=16)
    ax1.set_ylabel(gylabel,fontsize=16)
    plt.style.use('seaborn-colorblind')
    plt.tight_layout()

          
    # Plot n(r) data
    fig2, ax2 = plt.subplots()
    nxlabel = rf'$r$ (nm)'
    nylabel = rf'$n_{{\mathrm{{{rdf_refatom}-{rdfkey}}}}}(r)$'
    ax2.set_xlabel(nxlabel,fontsize=16)
    ax2.set_ylabel(nylabel,fontsize=16)
    plt.style.use('seaborn-colorblind')
    plt.tight_layout()

    # Plot combined g(r)/n(r) data
    fig3, ax3 = plt.subplots()
    gnxlabel = rf'$r$ (nm)'
    gnylabel = rf'$g_{{\mathrm{{{rdf_refatom}-{rdfkey}}}}}(r)$, $n_{{\mathrm{{{rdf_refatom}-{rdfkey}}}}}(r)$'
    ax3.set_xlabel(gnxlabel,fontsize=16)
    ax3.set_ylabel(gnylabel,fontsize=16)
    plt.style.use('seaborn-colorblind')
    plt.tight_layout()
    

    for inum,molval in enumerate(molality): # loop in molval

        print(f"Analyzing {molval} m")
        nfylcnt = 0

        anadir = resultdir + '/rdf_mol_' + str(molval)
        if not os.path.isdir(anadir):
            warnings.warn(f'{anadir} does not exist')
            continue

        # Check file
        gfname  = anadir + '/rdf_' + rdf_refatom + '_all.xvg'
        if not os.path.exists(gfname):
            warnings.warn(f"{gfname} does not exist in {anadir} ")
            continue

        # Check file
        nfname  = anadir + '/nrdf_' + rdf_refatom + '_all.xvg'
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

