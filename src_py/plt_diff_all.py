# Plots diffusivity for different concentration
# Fits VFT 
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import warnings
import aux_plt as aux

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

# Read the Excel file
dirkey     = 'diffusivity'
dirsuffix  = 'results_all_trial'+str(trialnum)
anadir     = '../../FeTFSI/analyzed_results/' +  dirsuffix + '/density_all'
figdir     = '../../FeTFSI/figures/' + dirsuffix +'/fetfsi_'+str(fecharge) +'/' + dirkey + '_all'

inpfile = f'{anadir}/density.xlsx'

if not os.path.isdir(anadir):
    raise RuntimeError(f'{anadir} not found')

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
y_f = df["d_F(cm2/s)"]
err_f = df["ed_F(cm2/s)"]

if not os.path.isdir(figdir):
    os.mkdir(figdir)

#-------------------------------------------------------------------------------
# Create scatter plot with error bars
fig1,ax1 = plt.subplots()
ax1.set_xlabel(r'Concentration (m)')
ax1.set_ylabel(r'Diffusivity (m$^2$/s)')
plt.style.use('seaborn-colorblind')
plt.tight_layout()

plt.errorbar(x, 0.0001*y_fe, yerr=0.0001*err_fe,color=clr_arr[0],\
             marker='o',markersize=8,capsize=4, linestyle='none',label='Fe')

plt.errorbar(x, 0.0001*y_f, yerr=0.0001*err_f,color=clr_arr[1],\
             marker='s',markersize=8,capsize=4, linestyle='none',label='F (TFSI)')

plt.legend()
fig1.savefig(figdir + '/diff_all_Fe' + str(fecharge) +
             '_trial_' + str(trialnum) + '.png',dpi=fig1.dpi)
#-------------------------------------------------------------------------------

# Create scatter plot with error bars in semilog scale
fig2,ax2 = plt.subplots()
ax2.set_xlabel(r'Concentration (m)')
ax2.set_ylabel(r'Diffusivity (m$^2$/s)')
plt.style.use('seaborn-colorblind')
plt.tight_layout()

plt.errorbar(x, 0.0001*y_fe, yerr=0.0001*err_fe,marker='o',\
             capsize=4, linestyle='none',label='Fe')

plt.errorbar(x, 0.0001*y_f, yerr=0.0001*err_f,marker='s',\
             capsize=4, linestyle='none',label='F (TFSI)')

plt.yscale('log')
plt.legend()
fig2.savefig(figdir + '/logdiff_all_Fe' + str(fecharge) +
             '_trial_' + str(trialnum) + '.png',dpi=fig1.dpi)

#-------------------------------------------------------------------------------

# Fit with VFT data

# Keep only positive diffusivities

mask_fe = (y_fe > 0)
mask_f  = (y_f > 0)

x_fe = x[mask_fe]
x_f  = x[mask_f]

y_fe = y_fe[mask_fe]
err_fe = err_fe[mask_fe]

y_f = y_f[mask_f]
err_f = err_f[mask_f]

# Fit Fe in log-space
popt_fe, pcov_fe = aux.curve_fit(log_vft_conc, x_fe,
                                 logy_fe,p0=p0_fe,\
                                 sigma=siglog_fe,\
                                 absolute_sigma=True,\
                                 bounds=(lower_bounds, upper_bounds),\
                                 maxfev=20000)

logD0_fe, B_fe, c0_fe = popt_fe
perr_fe = np.sqrt(np.diag(pcov_fe))
logD0_fe_err, B_fe_err, c0_fe_err = perr_fe

# Fit F in log-space
popt_f, pcov_f = aux.curve_fit(log_vft_conc, x_f, logy_f, p0=p0_f,\
                               sigma=siglog_f,absolute_sigma=True,\
                               bounds=(lower_bounds, upper_bounds),
                               maxfev=20000)

logD0_f, B_f, c0_f = popt_f
perr_f = np.sqrt(np.diag(pcov_f))
logD0_f_err, B_f_err, c0_f_err = perr_f

# Convert logD0 back to D0
D0_fe = np.exp(logD0_fe)
D0_f  = np.exp(logD0_f)

# Propagate uncertainty approximately: sigma(D0) = D0 * sigma(logD0)
D0_fe_err = D0_fe * logD0_fe_err
D0_f_err  = D0_f  * logD0_f_err

# Smooth curves
xfit = np.linspace(np.min(x), np.max(x), 400)

logfit_fe = log_vft_conc(xfit, *popt_fe)
logfit_f  = log_vft_conc(xfit, *popt_f)

yfit_fe = np.exp(logfit_fe)
yfit_f  = np.exp(logfit_f)

# Plot
fig3,ax3 = plt.subplots()
ax1.set_xlabel(r'Concentration (m)')
ax1.set_ylabel(r'Diffusivity (m$^2$/s)')
plt.style.use('seaborn-colorblind')
plt.tight_layout()

plt.errorbar(x, 0.0001*y_fe, yerr=0.0001*err_fe,marker='o',\
             color=clr_arr[0],markersize=8,capsize=4, \
             linestyle='none',label=r'Fe$^{3+}$')

plt.errorbar(x, 0.0001*y_f, yerr=0.0001*err_f,marker='s',\
             color = clr_arr[1],markersize=8,capsize=4, \
             linestyle='none',label='F (TFSI)')

plt.plot(xfit, yfit_fe,color=clr_arr[0],linestyle='--',\
         label=r'VFT Fit - Fe$^{3+}$)')
plt.plot(xfit, yfit_f,color=clr_arr[1],linestyle='--',\
         label=r'VFT fit - F(TFSI))')

plt.yscale("log")
plt.xlabel("Concentration (m)")
plt.ylabel("Diffusivity (m$^2$/s)")
plt.legend()
plt.tight_layout()
fig3.savefig(figdir + '/VFTfit_diff_all_Fe' + str(fecharge) +
             '_trial_' + str(trialnum) + '.png',dpi=fig3.dpi)



# Write fitted parameters
with open(f'{anadir}/D_fitted.dat','w') as fid:
    fid.write("Fe log-space fit parameters:")
    fid.write(f"log(D0) = {logD0_fe:.6f} Â± {logD0_fe_err:.6f}")
    fid.write(f"D0      = {D0_fe:.6e} Â± {D0_fe_err:.6e}")
    fid.write(f"B       = {B_fe:.6e} Â± {B_fe_err:.6e}")
    fid.write(f"c0      = {c0_fe:.6f} Â± {c0_fe_err:.6f}")
    
    fid.write("\nF log-space fit parameters:")
    fid.write(f"log(D0) = {logD0_f:.6f} Â± {logD0_f_err:.6f}")
    fid.write(f"D0      = {D0_f:.6e} Â± {D0_f_err:.6e}")
    fid.write(f"B       = {B_f:.6e} Â± {B_f_err:.6e}")
    fid.write(f"c0      = {c0_f:.6f} Â± {c0_f_err:.6f}")
