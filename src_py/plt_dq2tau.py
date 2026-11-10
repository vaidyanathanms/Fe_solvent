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


markers = ['o', 's', '^', 'D', 'v', 'P', 'X', '*']

# Inputs
fecharge = 3
trialnum = 5
iontype  = ['ion']
molality = ['0.1','0.2','0.4','1.0','2.5','5.0']
only_plot = 0 # 1 if just to plot

# Directory details
dirkey     = 'tau'
dirsuffix  = 'results_all_trial'+str(trialnum)+'/fetfsi_'+str(fecharge)
headdir    = '../../FeTFSI/' + dirsuffix
figdir     = '../../FeTFSI/figures/' + dirsuffix +'/' + dirkey + '_all'
anadir_tau = '../../FeTFSI/analyzed_results/' + dirsuffix + '/tau_results'
anadir_diff = '../../FeTFSI/analyzed_results/' + dirsuffix + '/dens_diff_all'


# Diff-data
diff_inpfile    = f'{anadir_diff}/dens_and_diff.xlsx'

if not os.path.isdir(anadir_tau):
    raise RuntimeError(f'analysis directory: {anadir_tau} not found')
if not os.path.isdir(anadir_diff):
    raise RuntimeError(f'analysis directory: {anadir_diff} not found')
if not os.path.exists(diff_inpfile):
    raise RuntimeError(f'{diff_inpfile} not found in {anadir_diff}')


# Output file
foutname = anadir_tau + '/diff_tau_q_calc_new'
col_names = ['conc_m','dval_cm2sinv','qval_nminv','beta','tauavg_ps',\
             'tauavginv_psinv','q^2D_psinv','q^2Dtau']
#--------------Main analysis------------------------------------------

if only_plot == 0:
    # Extract columns from diff_data
    df_diff = pd.read_excel(diff_inpfile,sheet_name='Fe'+str(fecharge)+
                            '_trial_'+str(trialnum))
    xmol = df_diff["mol"]
    diff_fe = df_diff["d_Fe(cm2/s)"]
    diff_tfsi = df_diff["d_TFSI_COM(cm^2/s)"]

    # Load tau data
    for atype in iontype:

        # Initialize output arrays
        q_arr      = np.array([])
        beta_arr   = np.array([])
        tauavg_arr = np.array([])
        diff_arr   = np.array([])
        conc_arr   = np.array([])
        tauinv_arr = np.array([])
        q2D_arr    = np.array([])
        q2dtau_arr = np.array([])
        
        for inum,molval in enumerate(molality): # loop in runarr
            
            Fe_diffval = diff_fe.loc[xmol == float(molval)].item()

            print(f"Analyzing {atype} tau for {molval}")
            
            if float(molval) > 1.0:
                rtype = 'main'
            else:
                rtype = 'smallrun'
            
            # Check file
            list_fnames = glob.glob(f'{anadir_tau}/tau{atype}*{rtype}*{molval}*')
            if list_fnames == []:
                warnings.warn(f'No tau{atype}*{rtype}* files exist in {anadir}')
                continue
            fname = max(list_fnames,key=os.path.getmtime)

            data = np.loadtxt(fname,comments = "#")
            
            qval = data[:,0]
            beta = data[:,2]
            tau_avg = data[:,3]

            mask = np.isfinite(qval) & np.isfinite(tau_avg)
        
            qval = qval[mask]
            beta = beta[mask]
            tau_avg = tau_avg[mask]

            locdiff_arr = 100*np.full(len(qval),Fe_diffval) # convert to nm^2/ps
            locconc_arr = np.full(len(qval),float(molval))
        
            tauavg_inv = 1/tau_avg
            qsquaredD  = qval*qval*locdiff_arr
            qsqDtau    = qsquaredD*tau_avg

            q_arr      = np.append(q_arr,qval)
            beta_arr   = np.append(beta_arr,beta)
            tauavg_arr = np.append(tauavg_arr,tau_avg)
            diff_arr   = np.append(diff_arr,locdiff_arr)
            conc_arr   = np.append(conc_arr,locconc_arr)
            tauinv_arr = np.append(tauinv_arr,tauavg_inv)
            q2D_arr    = np.append(q2D_arr,qsquaredD)
            q2dtau_arr = np.append(q2dtau_arr,qsqDtau)
            
            
    combined_data = np.column_stack((conc_arr,diff_arr,q_arr,beta_arr,\
                                     tauavg_arr,tauinv_arr,q2D_arr,\
                                     q2dtau_arr))

    df_out = pd.DataFrame(combined_data, columns=col_names)
    df_out.to_excel(f"{foutname}.xlsx", index=False)


#------------------------Plot outputs---------------------------------------
fname = anadir_tau + '/diff_tau_q_calc_new.xlsx'
if not os.path.exists(fname):
    raise RuntimeError(f'{fname} does not exist')


# Plot q^2D vs 1/tau
fig1,ax1 = plt.subplots()
df = pd.read_excel(fname,usecols=['conc_m','q^2D_psinv','tauavginv_psinv'])
avgdf = df.groupby(['conc_m','q^2D_psinv'],as_index = False).\
    agg(tauinvmean=('tauavginv_psinv','mean'),\
        tauinvstd=('tauavginv_psinv','std'),\
        n=('tauavginv_psinv','count'))

fitx = np.linspace(df['q^2D_psinv'].min(),df['q^2D_psinv'].max(),50)
trendline_coeff = np.polyfit(df['q^2D_psinv'].to_numpy(),\
                             df['tauavginv_psinv'].to_numpy(),1)
trendline_function = np.poly1d(trendline_coeff)

for i, (conc, group) in enumerate(avgdf.groupby('conc_m', sort=True)):
    points = ax1.scatter(group['q^2D_psinv'],group['tauinvmean'],\
                        marker=markers[i % len(markers)],
                        s=100,color=clr_arr[i],label=f'{conc} m')
    repeated = group.loc[group["n"]>1]
    if not repeated.empty:
        ax1.errorbar(group['q^2D_psinv'],group['tauinvmean'],
            yerr=repeated["tauinvstd"],fmt="none",
            ecolor=points.get_facecolors()[0],
            capsize=4,
        )

plt.plot(fitx, trendline_function(fitx), color="black", ls="--",lw=2)        
plt.style.use('seaborn-colorblind')
plt.legend()
plt.tight_layout()
ax1.set_xlabel(r'$q^{2} D_{Fe^{3+}}$ (ps$^{-1}$)', fontsize=16)
ax1.set_ylabel(r'$1/\tau$ (ps$^{-1}$)', fontsize=16)
fig1.savefig(f'{figdir}/tauinvdiffq2.png',dpi=fig1.dpi)

# Plot q^D\tau vs q
fig2,ax2 = plt.subplots()
df = pd.read_excel(fname,usecols=['conc_m','qval_nminv','q^2Dtau'])
avgdf = df.groupby(['conc_m','qval_nminv'],as_index = False).\
    agg(q2dtmean=('q^2Dtau','mean'),\
        q2dtstd=('q^2Dtau','std'),\
        n=('q^2Dtau','count'))

fitx = np.linspace(df['qval_nminv'].min(),df['qval_nminv'].max(),50)
fity = np.full(len(fitx),1)
for i, (conc, group) in enumerate(avgdf.groupby('conc_m', sort=True)):
    points = ax2.scatter(group['qval_nminv'],group['q2dtmean'],\
                        marker=markers[i % len(markers)],
                        s=100,color=clr_arr[i],label=f'{conc} m')
    repeated = group.loc[group["n"]>1]
    if not repeated.empty:
        ax2.errorbar(group['qval_nminv'],group['q2dtmean'],
            yerr=repeated["q2dtstd"],fmt="none",
            ecolor=points.get_facecolors()[0],
            capsize=4,
        )

        
plt.plot(fitx,fity, color="black", ls="--",lw=2)
plt.style.use('seaborn-colorblind')
ax2.set_xlim([2,8])
#ax2.set_ylim([0.5,1.2])
ax2.set_ylabel(r'$q^{2} D_{Fe^{3+}} \tau$', fontsize=16)
ax2.set_xlabel(r'$q$ (nm$^{-1}$)', fontsize=16)
#plt.legend(loc='upper right')
plt.tight_layout()
fig2.savefig(f'{figdir}/productall.png',dpi=fig2.dpi)


