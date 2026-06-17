# Plot all diff curves from fit_vft_and_plt_diff.py output
# Use ONLY IF o/ps from fit_vft_and_plt_diff.py are present
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Color/line data; figure defaults
orange = '#FFA500'; dark_g = '#006400'; brown = '#8B4513'
clr_arr = [dark_g,orange,brown,'b','k','m']
mrk_arr = ['o','d','s']
lne_arr = ['-','--']
plt.rc('legend',fontsize=16) # fontsize of the legends
plt.rcParams.update({'font.size': 14}) # other fonts
plt.rcParams.update({'figure.autolayout': True})
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']

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

mask_fe = (y_fe > 0)
mask_tfsi  = (y_tfsi > 0)

x_fe = x[mask_fe]
x_tfsi  = x[mask_tfsi]

y_fe = y_fe[mask_fe]
err_fe = err_fe[mask_fe]

y_tfsi = y_tfsi[mask_tfsi]
err_tfsi = err_tfsi[mask_tfsi]


if not os.path.isdir(figdir):
    os.mkdir(figdir)
#-------------------------------------------------------------------------------
print("Plotting Fe/TFSI diffusivities with VFT fits...")
fig1,ax1 = plt.subplots()
ax1.set_xlabel(r'Concentration (m)',fontsize=16)
ax1.set_ylabel(r'Diffusivity (m$^2$/s)',fontsize=16)
plt.style.use('seaborn-colorblind')
plt.tight_layout()

xfit = np.linspace(np.min(x_fe), np.max(x_fe), 400)

logfit_fe = aux.log_vft_conc(xfit, *popt_fe)
logfit_f  = aux.log_vft_conc(xfit, *popt_f)

yfit_fe = np.exp(logfit_fe)
yfit_f  = np.exp(logfit_f)

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
ax1.set_xlabel("Concentration (m)",fontsize=16)
ax1.set_ylabel("Diffusivity (cm$^2$/s)",fontsize=16)
ax1.legend()
plt.tight_layout()
fig3.savefig(figdir + '/VFTfit_diff_all_Fe' + str(fecharge) +
             '_trial_' + str(trialnum) + '.png',dpi=fig3.dpi)





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
