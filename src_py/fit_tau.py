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
from scipy.optimize import curve_fit
import math

system = 'cades' # cades or lap

# Color/line data; figure defaults
sns.set_palette("colorblind")
orange = '#FFA500'; dark_g = '#006400'; brown = '#8B4513'
clr_arr = [dark_g,orange,brown,'b','k','m']
mrk_arr = ['o','d','s']
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
runtype  = ['smallrun','main']
iontype  = ['ion','countion']

# Directory/file details
if system == 'cades':
    dirsuffix  = '/lustre/or-scratch/cades-birthright/vm5/fetfsi/fetfsi_3/trial_5'
    resultdir  =  dirsuffix + '/results_all'
    figdir     =  dirsuffix + '/results_all/figtau_all'
    outheaddir =  dirsuffix + '/results_all/tau_results'
else:
    # Directory details
    dirkey     = 'tau'
    dirsuffix  = 'results_all_trial'+str(trialnum)+'/fetfsi_'+str(fecharge)
    headdir    = '../../FeTFSI/' + dirsuffix
    resultdir  = headdir + '/' + dirkey + '_all'
    figdir     = '../../FeTFSI/figures/' + dirsuffix +'/' + dirkey + '_all'
    outheaddir = '../../FeTFSI/analyzed_results/results_all_trial' + \
        str(trialnum) + 'fetfsi_' + str(fecharge) + '/tau_results'


if not os.path.isdir(figdir):
    os.mkdir(figdir)

if not os.path.isdir(outheaddir):
    os.mkdir(outheaddir)

if not os.path.isdir(resultdir):
    raise RuntimeError(f'{resultdir} not found')    

# Load data
for inum,molval in enumerate(molality): # loop in runarr

    if not system == 'cades':
        outdir   = '../../FeTFSI/analyzed_results/mol_' + str(molval)
    else:
        outdir   = outheaddir
        
    if not os.path.isdir(outdir):
        os.mkdir(outdir)

    print(f"Analyzing {molval}")
    
    anadir = resultdir + '/ana_mol_' + molval

    if not os.path.isdir(anadir):
        warnings.warn(f'{anadir} does not exist')
        continue

    for rtype in runtype:

        print(f"Analyzing {rtype} for {molval}")
        
        for atype in iontype:

            print(f"Analyzing {atype} for {molval}")
            
            # Check file
            list_fnames = glob.glob(f'{anadir}/Fskt{atype}*{rtype}*')
            if list_fnames == []:
                warnings.warn(f'No Fskt{atype}*{rtype}* files exist in {anadir}')
                continue

            # Initialize arrays
            q_arr    = np.array([])
            tau_arr  = np.array([]); tauavg_arr = np.array([])
            beta_arr = np.array([]); 
            tauerr_arr  = np.array([])
            betaerr_arr = np.array([]); 

            # Set output file and output figure data
            fig1,ax1 = plt.subplots()
            ax1.set_xlabel(r'$q$ (nm$^{-1}$)')
            ax1.set_ylabel(r'$F$_s$(q,t)$')
            #plt.style.use('seaborn-colorblind')
            plt.tight_layout()
    
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
                    q_arr    = np.append(q_arr,qval)
                    tau_arr  = np.append(tau_arr,np.nan)
                    beta_arr = np.append(beta_arr,np.nan)
                    tauerr_arr  = np.append(tauerr_arr,np.nan)
                    betaerr_arr = np.append(betaerr_arr,np.nan)
                    continue

                # Fit initial guesses
                tau0 = time_data[len(time_data)//5]
                beta0 = 0.8
                init_guess = [tau0,beta0]
                bounds = ([1e-10, 0.01],
                          [np.inf, 1.0]) # same order as init_guess

                # Fit
                popt, pcov = curve_fit(aux.stretched_exp,time_data,fskt_data,\
                                           p0=init_guess,bounds=bounds,maxfev=100000)
                tau_fit, beta_fit,  = popt
                tau_avg = tau_fit * math.gamma(1.0 + 1.0/beta_fit)
                perr = np.sqrt(np.diag(pcov))
                
                # Save to arrays
                q_arr      = np.append(q_arr,qval)
                tau_arr    = np.append(tau_arr,tau_fit)
                beta_arr   = np.append(beta_arr,beta_fit)
                tauavg_arr = np.append(tauavg_arr,tau_avg)
                tauerr_arr  = np.append(tauerr_arr,perr[0])
                betaerr_arr = np.append(betaerr_arr,perr[1])
                
                # Plot data
                ax1.plot(time_data,fskt_data,"o",markersize=2,\
                         color=clr_arr[inum],label='_data')
                tfit = np.linspace(time_data.min(),time_data.max(), 1000)
                yfit = aux.stretched_exp(tfit, *popt)
                ax1.plot(tfit,yfit,"-",linewidth=2,\
                         color=clr_arr[inum],label=str(molval) + ' m')

        
            # Save output file
            combined_data = np.column_stack((q_arr,tau_arr,beta_arr,\
                                             tauerr_arr,betaerr_arr))
            
            np.savetxt(f'{outdir}/tau{atype}_{rtype}_mol_{molval}.dat',\
                       combined_data,delimiter='\t',\
                       header='q\ttau\tbeta\ttau_avg\ttau_err\tbeta_err\n')
                
            plt.legend(loc=0)
            plt.tight_layout()
            fig1.savefig(f'{figdir}/tau{atype}_{rtype}_mol_{molval}.png',\
                         dpi=fig1.dpi)
            plt.close(fig1)
