# Plots g(r) and n(r)
# Import modules
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import warnings
import aux_plt as aux
import glob
import warnings
from pathlib import Path


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
fecharge = 3
trialnum = 5
molality = ['0.1','0.2','0.4','1.0','2.5','5.0']

# Plot keys
sirdf    = 0
simsd    = 0
sifsqt   = 0
sicaged  = 1

if sirdf:
    # RDF Nitrogen data


    # RDF inputs
    rdf_refatom = 'Fe'
    rdf_keys    = {'Fe':1,'NTFSI':2}
    maxxlimplt  = [3.0,1.2,1.2,1.0,0.3]
    maxylimplt = [5,5,15,15,10]
    rcol  = 0

    # Directory info
    dirkey     = 'rdf'
    dirsuffix  = 'results_all_trial'+str(trialnum)+'/fetfsi_'+str(fecharge)
    headdir    = '../../FeTFSI/' + dirsuffix
    resultdir  = headdir + '/' + dirkey + '_all'
    figdir     = '../../FeTFSI/figures/' + dirsuffix +'/' + 'si_all'

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


if simsd:

    # Directory/file details
    dirkey     = 'diff'
    dirsuffix  = 'results_all_trial'+str(trialnum) + '/fetfsi_' + str(fecharge)
    anadir     = '../../FeTFSI/' +  dirsuffix + '/' + dirkey + '_all'
    figdir     = '../../FeTFSI/figures/' + dirsuffix  +'/' + 'si_all'
    
    
    if not os.path.isdir(anadir):
        raise RuntimeError(f'analysis directory: {anadir} not found')

    # Plot msd data
    fig1, ax1 = plt.subplots()
    mxlabel = r'Time (ps)'
    mylabel = r'MSD ($\AA^{2}$)'
    ax1.set_xlabel(mxlabel,fontsize=16)
    ax1.set_ylabel(mylabel,fontsize=16)
    plt.style.use('seaborn-colorblind')
    plt.tight_layout()

    
    for inum,molval in enumerate(molality): # loop in molval
        
        print(f"Analyzing MSD of {molval} m")
        nfylcnt = 0

        # Check file - Li
        liname  = anadir + '/diff_' + molval + '.xvg'
        if not os.path.exists(liname):
            warnings.warn(f"{liname} does not exist in {anadir} ")
            continue
         
        # Check file - TFSI COM
        tfsiname  = anadir + '/diff_TFSI_com_' + molval + '.xvg'
        if not os.path.exists(tfsiname):
            warnings.warn(f"{tfsiname} does not exist in {anadir} ")
            continue

        
        # Open and parse Li-MSD file
        lialldata = aux.read_xvg_filedata(liname)
        # Open and parse TFSI-MSD file
        tfsialldata = aux.read_xvg_filedata(tfsiname)

        tlendata = int(0.7*len(lialldata[:,0]))
        tlidata   = lialldata[0:tlendata,0]
        limsddata = lialldata[0:tlendata,1]
        ttfsidata = tfsialldata[0:tlendata,0]
        tfsimsddata = tfsialldata[0:tlendata,1]
        
        ax1.plot(tlidata,100*limsddata,color=clr_arr[inum],linestyle='-',label=molval + ' m')
        ax1.plot(ttfsidata,100*tfsimsddata,color=clr_arr[inum],linestyle='--',label='_m')

    ax1.legend()
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    fig1.savefig(figdir + '/msd_all.png',dpi=fig1.dpi)
    plt.close(fig1)
    


if sifsqt:

    iontype  = ['ion','countion']

    
    # Directory/file details
    dirkey     = 'neigh_diff'
    dirsuffix  = 'results_all_trial'+str(trialnum)+'/fetfsi_'+str(fecharge)
    headdir    = '../../FeTFSI/' + dirsuffix
    resultdir  = headdir + '/' + dirkey + '_all'
    figdir     = '../../FeTFSI/figures/' + dirsuffix +'/' + 'si_all'


    if not os.path.isdir(figdir):
        os.mkdir(figdir)

    if not os.path.isdir(resultdir):
        raise RuntimeError(f'{resultdir} not found')    


    # Load data
    for inum,molval in enumerate(molality): # loop in runarr

        print(f"Analyzing {molval}")
        
        anadir = resultdir + '/ana_mol_' + molval

        # Set output file and output figure data
        fig1,ax1 = plt.subplots()
        ax1.set_xlabel(r'Time (ps)')
        ax1.set_ylabel(r'$F_s$ ($q,t$)')
        plt.style.use('seaborn-colorblind')
        plt.tight_layout()

        for atype in iontype:

            if float(molval) > 1.0:
                rtype = 'main'
            else:
                rtype = 'smallrun'
                
            print(f"Analyzing {atype} for {molval}")
            
            # Check file
            list_fnames = glob.glob(f'{anadir}/Fskt{atype}*{rtype}*')
            if list_fnames == []:
                warnings.warn(f'No Fskt{atype}*{rtype}* files exist in {anadir}')
                continue

        
            # Loop through the files
            for fname in list_fnames:
                # Get qval
                qval = float(Path(fname).name.split('_')[1])
                # Open and parse file
                with open(fname) as fin:
                    lines = (line.lstrip() for line in fin \
                             if not line.lstrip().startswith('#') and \
                             not line.lstrip().startswith('@'))
                    data  = np.loadtxt(lines,skiprows=1)
                    tdata = data[:,0]
                    fskt  = data[:,1]
        
                # Smooth to avoid repeated values
                time_smooth, fskt_smooth = aux.compress_repeated_y(tdata,fskt,\
                                                                       tolerance=0.0,\
                                                                       method="middle")

        
                fskt_smooth /= fskt_smooth[0]
                
                # Truncate to the first negative value
                time_data, fskt_data = aux.truncate(time_smooth,fskt_smooth)
                
                
                # Check if y-values are going higher
                if len(fskt_data[fskt_data > 1.0]) != 0:
                    warnings.warn(f'Fskt values are not monotonically decreasing for {qval}')
                        
                # Plot data
                ax1.plot(time_data,fskt_data,"o",markersize=2,\
                         color=clr_arr[inum],label='_data')
                tfit = np.linspace(time_data.min(),time_data.max(), 1000)
                yfit = aux.stretched_exp(tfit, *popt)
                ax1.plot(tfit,yfit,"-",linewidth=2,\
                         color=clr_arr[inum],label=r'$q = $ '+ str(qval)+ ' $\AA^{-1}$')

            
        ax1.legend()
        ax1.set_xscale('log')
        ax1.set_yscale('log')
        fig1.savefig(f'{figdir}/fsqt{atype}_{rtype}_mol_{molval}.png',dpi=fig1.dpi)
        plt.close(fig1)

    
if sicaged:

    cagetype  = ['caged','uncaged']

    
    # Directory/file details
    dirkey     = 'neigh_diff'
    dirsuffix  = 'results_all_trial'+str(trialnum)+'/fetfsi_'+str(fecharge)
    headdir    = '../../FeTFSI/' + dirsuffix
    resultdir  = headdir + '/' + dirkey + '_all'
    figdir     = '../../FeTFSI/figures/' + dirsuffix +'/' + 'si_all'


    if not os.path.isdir(figdir):
        os.mkdir(figdir)

    if not os.path.isdir(resultdir):
        raise RuntimeError(f'{resultdir} not found')    


    # Load data
    for inum,cagestr in enumerate(cagetype): # loop in runarr

        print(f"Analyzing {cagestr}")
        
        anadir = resultdir + '/fskt_' + cagestr

        # Set output file and output figure data
        fig1,ax1 = plt.subplots()
        ax1.set_xlabel(r'Time (ps)')
        ax1.set_ylabel(r'$F_s$ ($q,t$)')
        plt.tight_layout()
        curves = []
        files = glob.glob(anadir + '/Fsktion*')

        for ifyle,fyle in enumerate(files):
            path = Path(fyle)
            _,q_text,conc_text,_ = path.stem.split("_")
            qval = float(q_text); conc = float(conc_text)
            with open(fyle,'r') as fin:
                time,fs = np.loadtxt(fin,skiprows=1,unpack=True)
            curves.append((qval,conc,time,fs))

        curves.sort(key=lambda curve: (curve[1], curve[0]))

        colors = plt.get_cmap("tab10").colors

        for fid, (qval, conc, time, fs) in enumerate(curves):
            positive = time > 0
            ax1.semilogx(
                time[positive],
                fs[positive],
                marker=mrk_arr[fid],
                markersize=3,
                linestyle=None,
                linewidth=1.5,
                color=colors[fid % len(colors)],
                label=rf"$q = {qval:g}\ \mathrm{{nm}}^{{-1}}$; {conc:g} m",
            )

                
        ax1.legend(loc='lower left',fontsize=12)
        fig1.savefig(f'{figdir}/fsqt_{cagestr}.png',dpi=fig1.dpi)
        plt.close(fig1)
