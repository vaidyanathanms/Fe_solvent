import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
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
plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']

# Inputs
fecharge = 3
trialnum = 5

# Read the Excel file
inpfile = '../../FeTFSI/analyzed_results/density_all/density.xlsx'
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

figdir   = '../../FeTFSI/figures/diff_all'

if not os.path.isdir(figdir):
    os.mkdir(figdir)

# Keep only positive diffusivities
# ----------------------------
mask_fe = (y_fe > 0)
mask_f  = (y_f > 0)

x_fe = x[mask_fe]
x_f  = x[mask_f]

y_fe = y_fe[mask_fe]
err_fe = err_fe[mask_fe]

y_f = y_f[mask_f]
err_f = err_f[mask_f]

# ----------------------------
# VFT-like model in concentration
# D(c) = D0 * exp[-B / (c0 - c)]
# log(D) = log(D0) - B / (c0 - c)
# We fit log(D) directly
# ----------------------------
def log_vft_conc(c, log10D0, B, c0):
    return log10D0 - (B / (c0 - c)) / np.log(10)

# log-10 of diffusivity
logy_fe = np.log10(y_fe)
logy_f  = np.log10(y_f)

# Propagated uncertainty in log(D): sigma_logD ~ sigma_D / D
# Put a fallback if any error bar is zero or missing
eps = 1e-30
siglog_fe = np.where(err_fe > 0, err_fe / np.maximum(y_fe, eps), 1.0)
siglog_f  = np.where(err_f  > 0, err_f  / np.maximum(y_f,  eps), 1.0)

xmax = np.max(x)

# Initial guesses
p0_fe = [np.log10(np.max(y_fe)), 1.0, xmax + 1.0]
p0_f  = [np.log10(np.max(y_f)),  1.0, xmax + 1.0]

# Bounds: c0 must be > max(x)
lower_bounds = [-np.inf, 0.0, xmax + 1e-6]
upper_bounds = [ np.inf, np.inf, np.inf]

# ----------------------------
# Fit Fe in log-space
# ----------------------------
popt_fe, pcov_fe = curve_fit(
    log_vft_conc, x_fe, logy_fe,
    p0=p0_fe,
    sigma=siglog_fe,
    absolute_sigma=True,
    bounds=(lower_bounds, upper_bounds),
    maxfev=20000
)

logD0_fe, B_fe, c0_fe = popt_fe
perr_fe = np.sqrt(np.diag(pcov_fe))
logD0_fe_err, B_fe_err, c0_fe_err = perr_fe

# ----------------------------
# Fit F in log-space
# ----------------------------
popt_f, pcov_f = curve_fit(
    log_vft_conc, x_f, logy_f,
    p0=p0_f,
    sigma=siglog_f,
    absolute_sigma=True,
    bounds=(lower_bounds, upper_bounds),
    maxfev=20000
)

logD0_f, B_f, c0_f = popt_f
perr_f = np.sqrt(np.diag(pcov_f))
logD0_f_err, B_f_err, c0_f_err = perr_f

# Convert logD0 back to D0
D0_fe = np.exp(logD0_fe)
D0_f  = np.exp(logD0_f)

# Propagate uncertainty approximately: sigma(D0) = D0 * sigma(logD0)
D0_fe_err = D0_fe * logD0_fe_err
D0_f_err  = D0_f  * logD0_f_err

# ----------------------------
# Smooth curves
# ----------------------------
xfit = np.linspace(np.min(x), np.max(x), 400)

logfit_fe = log_vft_conc(xfit, *popt_fe)
logfit_f  = log_vft_conc(xfit, *popt_f)

yfit_fe = np.exp(logfit_fe)
yfit_f  = np.exp(logfit_f)

# ----------------------------
# Plot
# ----------------------------
plt.figure(figsize=(7, 5))


plt.errorbar(
    x_fe, y_fe, yerr=err_fe,
    fmt='o', linestyle='none', capsize=4,
    label='Fe data'
)
plt.errorbar(
    x_f, y_f, yerr=err_f,
    fmt='s', linestyle='none', capsize=4,
    label='F data'
)

plt.plot(
    xfit, yfit_fe,
    label=fr'Fe log-VFT fit ($c_0={c0_fe:.3f}\pm{c0_fe_err:.3f}$)'
)
plt.plot(
    xfit, yfit_f,
    label=fr'F log-VFT fit ($c_0={c0_f:.3f}\pm{c0_f_err:.3f}$)'
)

plt.yscale("log")
plt.xlabel("mol")
plt.ylabel("Diffusivity (m$^2$/s)")
plt.title("mol vs diffusivity with log-space VFT-like fits")
plt.legend()
plt.tight_layout()
plt.show()

# ----------------------------
# Print fitted parameters
# ----------------------------
print("Fe log-space fit parameters:")
print(f"log(D0) = {logD0_fe:.6f} Â± {logD0_fe_err:.6f}")
print(f"D0      = {D0_fe:.6e} Â± {D0_fe_err:.6e}")
print(f"B       = {B_fe:.6e} Â± {B_fe_err:.6e}")
print(f"c0      = {c0_fe:.6f} Â± {c0_fe_err:.6f}")

print("\nF log-space fit parameters:")
print(f"log(D0) = {logD0_f:.6f} Â± {logD0_f_err:.6f}")
print(f"D0      = {D0_f:.6e} Â± {D0_f_err:.6e}")
print(f"B       = {B_f:.6e} Â± {B_f_err:.6e}")
print(f"c0      = {c0_f:.6f} Â± {c0_f_err:.6f}")
