# Analysis and plot scripts for plotting
# 1. Anion neighbor distribution
# 2. CIP/SSIP and Free ions
# 3. Cluster size distribution
# Version: June-05-2026
# Author: Vaidyanathan M. Sethuraman

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import glob
import warnings
from pathlib import Path
import aux_plt as aux

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
molality = ['0.1','0.2','0.4','1.0','2.5','5.0']

# Flags
ionpair_flag = 0 # CIP/SSIP/free-plot
neigh_flag   = 0; neigh_keys  = {'Fe-An':3,'An-Fe':5}; maxneigh = 40 #neigh-plot
clust_flag   = 1; clust_keys  = {'Fe-An':2}; maxclust = 19500 #clust-plot

# File prefixes
ionpair_pref = 'conepairs_summary'
neigh_pref   = 'catanneigh_16_*main_nojump'
clust_pref   = 'clust_16_*main_nojump'

# Directory info
dirsuffix  = 'results_all_trial'+str(trialnum)+'/fetfsi_'+str(fecharge)
headdir    = '../../FeTFSI/' + dirsuffix
resultdir  = headdir + '/neigh_diff_all'
figdir     = '../../FeTFSI/figures/' + dirsuffix +'/neigh_all'
outheaddir = '../../FeTFSI/analyzed_results/' +  dirsuffix + '/neigh_results'

if not os.path.isdir(headdir):
    raise RuntimeError(f'No simulation head directory: {headdir} found')

if not os.path.isdir(figdir):
    os.mkdir(figdir)

if not os.path.isdir(outheaddir):
    os.mkdir(outheaddir)


# Load data
if ionpair_flag:
    
    cip_arr      = np.zeros(len(molality))
    ssip_arr     = np.zeros(len(molality))
    free_arr     = np.zeros(len(molality))
    mean_rad_arr = np.zeros(len(molality))
    for inum,molval in enumerate(molality): 

        print (f'Ionpair analysis for {molval} m')
    
        workdir = resultdir + '/ana_mol_' + molval
        
        if not os.path.isdir(workdir):
            warnings.warn(f'{workdir} does not exist')
            continue

        # Check file
        fname = workdir + '/' + ionpair_pref + '.csv'
        if not os.path.exists(fname):
            warnings.warn(f'No {fname} file exist in {workdir}')
            continue

        cip,ssip,free,meanr = aux.extract_ionpair_data(fname)

        cip_arr[inum] = cip; ssip_arr[inum] = ssip
        free_arr[inum] = free; mean_rad_arr[inum] = meanr
        
    aux.plot_grouped_fractions(molality,free_arr,ssip_arr,cip_arr,\
                               figfile=figdir + '/ionpair_dist.png')

        
#-------------------------------------------------------------------------------
if neigh_flag:
    
    for nid, (neighkey,colval) in enumerate(neigh_keys.items()):

        # Plot data
        fig, ax = plt.subplots()
        if nid == 0:
            aux.set_axes(ax,plt,r'Number of O-TFSI neighbors, $s$',r'$f_{\mathrm{neigh}}(s)$')
        else:
            aux.set_axes(ax,plt,r'Number of Fe$^{3+}$ neighbors, $s$',r'$f_{\mathrm{neigh}}(s)$')
        bar_width = 0.2
                
        all_data = np.zeros((maxneigh,len(molality)+1))

        for inum,molval in enumerate(molality): 

            workdir = resultdir + '/ana_mol_' + molval
        
            if not os.path.isdir(workdir):
                warnings.warn(f'{workdir} does not exist')
                continue

            print (f'Neighbor analysis for {neighkey} in {workdir}')

            xarr,yarr = aux.return_neigh_arrays(workdir,\
                                                neigh_pref,[0,colval-1],\
                                                0,maxneigh)
            if not len(xarr) > 1:
                continue
            
            if inum == 0:
                all_data[:,inum] = xarr - 1
            all_data[:,inum+1]   = yarr

        for i in range(len(molality)):
            ax.plot(all_data[:,0],all_data[:,i+1]/100,marker=mrk_arr[i],markersize=6,\
                    color=clr_arr[i],linestyle = '--',linewidth=2,\
                    label=f'{molality[i]} m')


        ax.set_xlim([-0.2, maxneigh+1])
        ax.set_ylim([0.0,0.35])
        ax.legend()
        fig.savefig(f'{figdir}/neigh_{neighkey}.png',dpi = fig.dpi)
        fig.savefig(f'{figdir}/neigh_{neighkey}.eps',format = 'eps')
        
#-------------------------------------------------------------------------------
if clust_flag:

    for nid, (clustkey,colval) in enumerate(clust_keys.items()):

        bar_width = 0.2
               
        # Plot data
        fig, ax = plt.subplots()
        aux.set_axes(ax,plt,r'Cluster size, $n$',r'$S$($n$)')
        
        fig2, ax2 = plt.subplots()
        aux.set_axes(ax,plt,r'Cluster size, $n$',r'$S$($n$)')
        
        for inum,molval in enumerate(molality): 

            workdir = resultdir + '/ana_mol_' + molval
        
            if not os.path.isdir(workdir):
                warnings.warn(f'{workdir} does not exist')
                continue

            print (f'Cluster analysis for {clustkey} in {workdir}')
            
            xarr,yarr = aux.return_clust_arrays(workdir,\
                                                 clust_pref,[0,colval-1],\
                                                 0,0)

            if not len(xarr) > 1:
                continue

            xarr = xarr[yarr>0.001]
            yarr = yarr[yarr>0.001]

            ax.plot(xarr,yarr,marker=mrk_arr[inum],markersize=6,linestyle='None',\
                    color=clr_arr[inum],label=f'{molality[inum]} m')
            ax2.plot(xarr,yarr,marker=mrk_arr[inum],markersize=6,linestyle='None',\
                    color=clr_arr[inum],label=f'{molality[inum]} m')
            
#        ax.set_xlim([0.99,100])
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.legend(loc='upper center')
        fig.savefig(f'{figdir}/clustall.png',dpi = fig.dpi)
        fig.savefig(f'{figdir}/clustall.eps',format = 'eps')
            
        ax2.set_xlim([0, 30])
        ax2.legend()
        fig2.savefig(f'{figdir}/clustinset.png',dpi = fig.dpi)
        fig2.savefig(f'{figdir}/clustinset.eps',format = 'eps')
            

        
