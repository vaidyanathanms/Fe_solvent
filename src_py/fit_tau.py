# To fit tau for autocorrelation functions
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import glob
import warnings
from pathlib import Path
import aux_plt as aux
from scipy.optimize import curve_fit

# Color/line data; figure defaults
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

# Load data
for inum,molval in enumerate(molality): # loop in runarr

    outdir   = '../../FeTFSI/analyzed_results/mol_' + str(molval)
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
        
    print(f"Analyzing {molval}")
    
    anadir = resultdir + '/mol_' + molval

    if not os.path.isdir(anadir):
        warnings.warn(f'{anadir} does not exist')
        continue
            
    # Check file
    list_fnames = glob.glob(anadir + '/Fsktion*small*')
    if list_fnames == []:
        warnings.warn(f'No Fskt small run files exist in {anadir}')
        continue

    # Initialize arrays
    q_arr    = np.array([])
    tau_arr  = np.array([]); tauavg_arr = np.array([])
    beta_arr = np.array([]); c0_arr = np.array([])
    tauerr_arr  = np.array([])
    betaerr_arr = np.array([]); c0err_arr = np.array([])

    # Set output file and output figure data
    fig1,ax1 = plt.subplots()
    ax1.set_xlabel(r'$q$ (nm$^{-1}$)')
    ax1.set_ylabel(r'$F$_s$(q,t)$')
    plt.style.use('seaborn-colorblind')
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
        
        # Check if y-values are going higher
        if len(fskt_smooth[fskt_smooth > 1.0]) != 0:
            warnings.warn(f'Fskt values are not monotonically decreasing for {qval}')
            q_arr    = np.append(q_arr,qval)
            tau_arr  = np.append(tau_arr,np.nan)
            beta_arr = np.append(beta_arr,np.nan)
            c0_arr   = np.append(c0_arr,np.nan)
            tauerr_arr  = np.append(tauerr_arr,np.nan)
            betaerr_arr = np.append(betaerr_arr,np.nan)
            c0err_arr   = np.append(c0err_arr,np.nan)
            continue

        # Truncate to the first negative value
        time_data, fskt_data = aux.truncate(time_smooth,fskt_smooth)

        # Fit initial guesses
        c0 = fskt_smooth[-1]
        tau0 = time_data[len(time_data)//5]
        beta0 = 0.8
        init_guess = [tau0,beta0,c0]
        bounds = ([1e-10, 0.01, 0.0],
                  [np.inf, 2.0, np.inf]) # same order as init_guess

        # Fit
        popt, pcov = aux.curve_fit(stretched_exp,time_data,fskt_data,\
                                   p0=p0,bounds=bounds,maxfev=100000)
        tau_fit, beta_fit, c0_fit = popt
        tau_avg = tau_fit * gamma(1.0 + 1.0/beta_fit)
        perr = np.sqrt(np.diag(pcov))

        # Save to arrays
        q_arr      = np.append(q_arr,qval)
        tau_arr    = np.append(tau_arr,tau_fit)
        beta_arr   = np.append(beta_arr,beta_fit)
        c0_arr     = np.append(c0_arr,c0_fit)
        tauavg_arr = np.append(c0_arr,tau_avg)
        tauerr_arr  = np.append(tauerr_arr,perr[0])
        betaerr_arr = np.append(betaerr_arr,perr[1])
        c0err_arr   = np.append(c0err_arr,perr[0])
                
        # Plot data
        ax1.plot(time_smooth,fskt_smooth,"o",markersize=2,\
                 color=clr_arr[inum],label='_data')
        tfit = np.linspace(time_smooth.min(),time_smooth.max(), 1000)
        yfit = aux.stretched_exp(tfit, *popt)
        ax1.plot(tfit,yfit,"-",linewidth=2,\
                 color=clr_arr[inum],label=str(molval) + ' m')

        
    # Save output file
    combined_data = np.column_stack((q_arr,tau_arr,beta_arr,c0_arr,\
                                     tauerr_arr,betaerr_arr,c0err_arr))
    
    np.savetxt(outdir + '/tausmall_mol_' + str(molval)+'.dat',\
               combined_data,delimiter='\t',\
               header='q\ttau\tbeta\tc0\ttau_avg\ttau_err\tbeta_err\tc0_err')

    plt.legend(loc=0)
    plt.tight_layout()
    fig1.savefig(figdir + '/tausmall_mol_' + str(molval)+ '.png',\
                 dpi=fig1.dpi)

