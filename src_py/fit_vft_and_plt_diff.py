# Plots diffusivity for different concentration
# Fits VFT 
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import warnings
import aux_plt as aux
import seaborn as sns

system   = 'cades' #cades or lap

# Color/line data; figure defaults
sns.set_palette("colorblind")
orange = '#FFA500'; dark_g = '#006400'; brown = '#8B4513'
clr_arr = [dark_g,orange,brown,'b','k','m']
mrk_arr = ['o','d','s']
lne_arr = ['-','--']
if system=='cades':
    font_dir = os.path.expanduser('~/.fonts')
plt.rc('legend',fontsize=16) # fontsize of the legends
plt.rcParams.update({'font.size': 14}) # other fonts
plt.rcParams.update({'figure.autolayout': True})
plt.rcParams['font.sans-serif'] = ['Liberation Sans', 'DejaVu Sans', 'Arial']
plt.rcParams['font.family'] = 'sans-serif'

# Inputs
fecharge = 3
trialnum = 5
dirkey   = 'diffusivity'

# Directory/file details
if system == 'cades':
    dirsuffix  = '/lustre/or-scratch/cades-birthright/vm5/fetfsi/fetfsi_3/trial_5/results_all'
    anadir     =  dirsuffix + '/dens_diff_all'
    figdir     =  dirsuffix + '/figs_all'
else:
    dirsuffix  = 'results_all_trial'+str(trialnum)
    anadir     = '../../FeTFSI/analyzed_results/' +  dirsuffix + '/density_all'
    figdir     = '../../FeTFSI/figures/' + dirsuffix +'/fetfsi_'+str(fecharge) +'/' + dirkey + '_all'
    
inpfile = f'{anadir}/dens_and_diff.xlsx'

if not os.path.isdir(anadir):
    raise RuntimeError(f'analysis directory: {anadir} not found')

if not os.path.exists(inpfile):
    raise RuntimeError(f'{inpfile} not found in {anadir}')

df = pd.read_excel(inpfile,sheet_name='Fe'+str(fecharge)+
                   '_trial_'+str(trialnum))

# print column names
print(df.columns)

# Extract columns
x = df["mol"]
y_fe = df["d_Fe(cm2/s)"]
err_fe = df["ed_Fe(cm2/s)"]
y_tfsi = df["d_TFSI_COM(cm^2/s)"]
err_tfsi = df["ed_TFSI_COM(cm^2/s)"]

if not os.path.isdir(figdir):
    os.mkdir(figdir)

#-------------------------------------------------------------------------------
# Create scatter plot with error bars
fig1,ax1 = plt.subplots()
ax1.set_xlabel(r'Concentration (m)',fontsize=16)
ax1.set_ylabel(r'Diffusivity (m$^2$/s)',fontsize=16)
plt.style.use('seaborn-colorblind')
plt.tight_layout()
print("Plotting Fe/TFSI diffusivities ...")
plt.errorbar(x, 0.0001*y_fe, yerr=0.0001*err_fe,color=clr_arr[0],\
             marker='o',markersize=8,capsize=4, linestyle='none',label=r'Fe$^{3+}$')

plt.errorbar(x, 0.0001*y_tfsi, yerr=0.0001*err_tfsi,color=clr_arr[1],\
             marker='s',markersize=8,capsize=4, linestyle='none',label=r'TFSI$^{-}$')

ax1.legend()
fig1.savefig(figdir + '/diff_all_Fe' + str(fecharge) +
             '_trial_' + str(trialnum) + '.png',dpi=fig1.dpi)
#-------------------------------------------------------------------------------

# Create scatter plot with error bars in semilog scale
fig2,ax2 = plt.subplots()
ax2.set_xlabel(r'Concentration (m)',fontsize=16)
ax2.set_ylabel(r'Diffusivity (m$^2$/s)',fontsize=16)
plt.style.use('seaborn-colorblind')
plt.tight_layout()
print("Plotting Fe/TFSI diffusivities in semilog scale...")
plt.errorbar(x, 0.0001*y_fe, yerr=0.0001*err_fe,marker='o',\
             capsize=4, linestyle='none',label=r'Fe$^{3+}$')

plt.errorbar(x, 0.0001*y_tfsi, yerr=0.0001*err_tfsi,marker='s',\
             capsize=4, linestyle='none',label=r'TFSI$^{-}$')

plt.yscale('log')
ax2.legend()
fig2.savefig(figdir + '/logdiff_all_Fe' + str(fecharge) +
             '_trial_' + str(trialnum) + '.png',dpi=fig1.dpi)

#-------------------------------------------------------------------------------

# Fit with VFT data

# Keep only positive diffusivities
print("Computing VFT fits ...")
mask_fe = (y_fe > 0)
mask_tfsi  = (y_tfsi > 0)

x_fe = x[mask_fe]
x_tfsi  = x[mask_tfsi]

y_fe = y_fe[mask_fe]
err_fe = err_fe[mask_fe]

y_tfsi = y_tfsi[mask_tfsi]
err_tfsi = err_tfsi[mask_tfsi]

# natural log of diffusivity
logy_fe    = np.log(y_fe)
logy_tfsi  = np.log(y_tfsi)

# Propagated uncertainty in log(D): sigma_logD ~ sigma_D / D
# Put a fallback if any error bar is zero or missing
eps = 1e-30
siglog_fe = np.where(err_fe > 0, err_fe / np.maximum(y_fe, eps), 1.0)
siglog_f  = np.where(err_tfsi  > 0, err_tfsi  / np.maximum(y_tfsi,  eps), 1.0)

xmax = np.max(x)

# Initial guesses
p0_fe = [np.log(np.max(y_fe)), 1.0, xmax + 1.0]
p0_f  = [np.log(np.max(y_tfsi)),  1.0, xmax + 1.0]

# Bounds: c0 must be > max(x)
lower_bounds = [-np.inf, 0.0, xmax + 1e-6]
upper_bounds = [ np.inf, np.inf, np.inf]

# Fit Fe in log-space
popt_fe, pcov_fe = aux.curve_fit(aux.log_vft_conc, x_fe,logy_fe,\
                                 p0=p0_fe,sigma=siglog_fe,\
                                 absolute_sigma=True,\
                                 bounds=(lower_bounds, upper_bounds),\
                                 maxfev=20000)

logD0_fe, B_fe, c0_fe = popt_fe
perr_fe = np.sqrt(np.diag(pcov_fe))
logD0_fe_err, B_fe_err, c0_fe_err = perr_fe

# Fit F in log-space
popt_f, pcov_f = aux.curve_fit(aux.log_vft_conc,x_tfsi,logy_tfsi,p0=p0_f,\
                               sigma=siglog_f,absolute_sigma=True,\
                               bounds=(lower_bounds, upper_bounds),
                               maxfev=20000)

logD0_f, B_f, c0_f = popt_f
perr_tfsi = np.sqrt(np.diag(pcov_f))
logD0_f_err, B_f_err, c0_f_err = perr_tfsi

# Convert logD0 back to D0
D0_fe = np.exp(logD0_fe)
D0_f  = np.exp(logD0_f)

# Propagate uncertainty approximately: sigma(D0) = D0 * sigma(logD0)
D0_fe_err = D0_fe * logD0_fe_err
D0_f_err  = D0_f  * logD0_f_err

# Smooth curves
xfit = np.linspace(np.min(x), np.max(x), 400)

logfit_fe = aux.log_vft_conc(xfit, *popt_fe)
logfit_f  = aux.log_vft_conc(xfit, *popt_f)

yfit_fe = np.exp(logfit_fe)
yfit_f  = np.exp(logfit_f)

# Plot
fig3,ax3 = plt.subplots()
plt.style.use('seaborn-colorblind')
plt.tight_layout()
print("Plotting Fe/TFSI diffusivities with VFT fits...")
plt.errorbar(x, y_fe, yerr=err_fe,marker='o',\
             color=clr_arr[0],markersize=8,capsize=4, \
             linestyle='none',label=r'Fe$^{3+}$')

plt.errorbar(x, y_tfsi, yerr=err_tfsi,marker='s',\
             color = clr_arr[1],markersize=8,capsize=4, \
             linestyle='none',label='TFSI$^{-}$')

plt.plot(xfit, yfit_fe,color=clr_arr[0],linestyle='--',\
         linewidth = 2.5, label=r'VFT Fit - Fe$^{3+}$')
plt.plot(xfit, yfit_f,color=clr_arr[1],linestyle='--',\
         linewidth = 2.5, label=r'VFT Fit - TFSI$^{-}$')

plt.yscale("log")
ax3.set_xlabel("Concentration (m)",fontsize=16)
ax3.set_ylabel("Diffusivity (cm$^2$/s)",fontsize=16)
ax3.legend()
plt.tight_layout()
fig3.savefig(figdir + '/VFTfit_diff_all_Fe' + str(fecharge) +
             '_trial_' + str(trialnum) + '.png',dpi=fig3.dpi)


print("Writing output to file..")
# Write fitted parameters
with open(f'{anadir}/D_fitted.dat','w') as fid:
    fid.write("Fe log-space fit parameters: \n")
    fid.write(f"log(D0) = {logD0_fe:.6f} ± {logD0_fe_err:.6f}\n")
    fid.write(f"D0      = {D0_fe:.6e} ± {D0_fe_err:.6e}\n")
    fid.write(f"B       = {B_fe:.6e} ± {B_fe_err:.6e}\n")
    fid.write(f"c0      = {c0_fe:.6f} ± {c0_fe_err:.6f}\n")
    
    fid.write("\nF log-space fit parameters:")
    fid.write(f"log(D0) = {logD0_f:.6f} ± {logD0_f_err:.6f}\n")
    fid.write(f"D0      = {D0_f:.6e} ± {D0_f_err:.6e}\n")
    fid.write(f"B       = {B_f:.6e} ± {B_f_err:.6e}\n")
    fid.write(f"c0      = {c0_f:.6f} ± {c0_f_err:.6f}\n")

# Plot D_Fe/D_F
fig4,ax4 = plt.subplots()
rat_data = df['d_TFSICOM/d_Fe']
ax4.set_xlabel(r'Concentration (m)',fontsize=16)
ax4.set_ylabel(r'$D_{\rm{TFSI}^{-}}$/$D_{\rm{Fe}^{3+}}$',fontsize=16)
plt.style.use('seaborn-colorblind')
plt.tight_layout()
print("Plotting ratio of Fe/TFSI diffusivities")
plt.plot(x, rat_data, marker='o',color=clr_arr[0],\
         markersize=8,linestyle = '--',linewidth=2.5)
fig4.savefig(figdir + '/rat_diff' + str(fecharge) +
             '_trial_' + str(trialnum) + '.png',dpi=fig3.dpi)
