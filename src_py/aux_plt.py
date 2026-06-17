# Auxiliary functions for analysis and plotting
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import csv
import glob
import warnings
from pathlib import Path
import re

# General copy script
def gencpy(dum_maindir,dum_destdir,fylname):

    srcfyl = dum_maindir + '/' + fylname

    if not os.path.exists(srcfyl):
        print('ERROR: ', srcfyl, 'not found')
        return

    desfyl = dum_destdir + '/' + fylname
    shutil.copy2(srcfyl,desfyl)
#-------------------------------------------------------------------------------

# Change width of seaborn barplot
def change_width(ax, new_value):
    for patch in ax.patches:
        current_width = patch.get_width()
        diff = current_width - new_value

        # Change the bar width
        patch.set_width(new_value)

        # Recenter the bar
        patch.set_x(patch.get_x() + diff * .5)
#-------------------------------------------------------------------------------

# Set axes labels
def set_axes(axhdl,plt,xlabel,ylabel):
    plt.style.use('seaborn-colorblind')
    plt.tight_layout()
    axhdl.set_xlabel(xlabel)
    axhdl.set_ylabel(ylabel)
    change_width(axhdl,0.2)
#-------------------------------------------------------------------------------

# find axes limits
def axlims(yminref,ymin,ymaxref,ymax):
    if ymax > ymaxref: ymaxref = ymax
    if ymin < yminref: yminref = ymin
    return yminref, ymaxref
#-------------------------------------------------------------------------------


# Fskt shows plateau at different times.
# But use just one representative time point per plateau
def compress_repeated_y(x, y, tolerance=0.0, method="middle"):
    x_new = []
    y_new = []

    start = 0
    n = len(y)

    for i in range(1, n + 1):
        #Check if adjacent y values are within a tolerance
        if i == n or abs(y[i] - y[start]) > tolerance:
            x_block = x[start:i]
            y_block = y[start:i]

            if method == "middle":
                idx = len(x_block) // 2
                x_rep = x_block[idx]
            elif method == "mean":
                x_rep = np.mean(x_block)
            elif method == "first":
                x_rep = x_block[0]
            elif method == "last":
                x_rep = x_block[-1]
            else:
                raise ValueError("method must be 'middle', 'mean', 'first', or 'last'")

            y_rep = np.mean(y_block)

            x_new.append(x_rep)
            y_new.append(y_rep)

            start = i

    return np.array(x_new), np.array(y_new)
#-------------------------------------------------------------------------------

# Stretched exponential function
# Fs(k,t)/Fs(k,0) = exp[-(t/tau)^beta] + C
# Use only normalized values
def stretched_exp(t, tau, beta, C):
    return np.exp(- (t / tau)**beta) + C
#-------------------------------------------------------------------------------

# Truncate data at the first point where Fskt goes below 0
# Keep only data before that point
def truncate(time,fskt):
    negative_indices = np.where(fskt < 0.0)[0]

    if len(negative_indices) > 0:
        first_negative = negative_indices[0]

        # keep data only before Fskt becomes negative
        time = time[:first_negative]
        fskt = fskt[:first_negative]

    return time,fskt
#-------------------------------------------------------------------------------

# R squared
def comp_rsquared(xdata, ydata, popt):
    residuals = ydata - stretched_exp(xdata, *popt)
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((ydata - np.mean(ydata))**2)
    r2 = 1 - ss_res / ss_tot
    return r2
#-------------------------------------------------------------------------------

# Extract data from conepairs_summary output
def extract_ionpair_data(inpfile):
    free_frac = None
    ssip_frac = None
    cip_frac  = None
    tfsi_effective_radius = None

    with open(inpfile, "r", newline="") as f:
        reader = csv.DictReader(f)

        for row in reader:
            metric = row["metric"].strip()

            if metric == "tfsi_effective_radius_A":
                tfsi_effective_radius = float(row["free"])

            elif metric == "raw_fractions":
                free_frac = float(row["free"])
                ssip_frac = float(row["SSIP"])
                cip_frac = float(row["CIP"])

    return cip_frac,ssip_frac,free_frac,tfsi_effective_radius
#-------------------------------------------------------------------------------

def plot_grouped_fractions(concentrations,free_arr,ssip_arr,cip_arr,\
                           figfile='../../figures/ssip.png'):
    x = np.arange(len(concentrations))   # group positions
    width = 0.25                         # width of each bar
    
    fig, ax = plt.subplots()
    
    ax.bar(x - width, free_arr, width, label="Free-ions")
    ax.bar(x, ssip_arr, width, label="SSIP")
    ax.bar(x + width, cip_arr, width, label="CIP")

    ax.set_xlabel("Concentration (m)")
    ax.set_ylabel("Fractions")
    ax.set_xticks(x)
    ax.set_xticklabels(concentrations)
    ax.set_ylim(0, 1.0)

    ax.legend()
    fig.tight_layout()

    plt.savefig(figfile, dpi=300)
    plt.close()
#-------------------------------------------------------------------------------

# Neighbor analysis
def return_neigh_arrays(simdir,strpref,colarr,nheaders=0,maxneigh=5):
    if len(colarr) != 2:
        raise RuntimeError("Expecting 2 columns for plotting neighbors")

    neigh_file = find_latest_file(glob.glob(simdir + '/' + strpref + '*'))
    if neigh_file == -1:
        warnings.warn(f'No neighbor file exists in {workdir}')
        return -1, -1

    sarr,neigharr = extract_data(neigh_file,colarr,nheaders,maxneigh)
    return sarr,neigharr
#-------------------------------------------------------------------------------

# Cluster analysis
def return_clust_arrays(simdir,strpref,colarr,nheaders=0,maxclust=20):
    if len(colarr) != 2:
        raise RuntimeError("Expecting 2 columns for plotting clusters")

    neigh_file = find_latest_file(glob.glob(simdir + '/' + strpref + '*'))
    if neigh_file == -1:
        warnings.warn(f'No cluster file exists in {workdir}')

    sarr,neigharr = extract_data(neigh_file,colarr,nheaders,maxclust)
    return sarr,neigharr
#-------------------------------------------------------------------------------

# Extract data from files according to columns
def extract_data(filename,colarr,skipl=1,refnr=0):
    data = np.loadtxt(filename, skiprows=skipl)  # Skips n header lines
    nrows = data.shape[0]
    if refnr == 0 or nrows < refnr:
        col1 = data[:, colarr[0]]  # Column 1
        col2 = data[:, colarr[1]]  # Column n
    else:
        col1 = data[:refnr, colarr[0]]  # Column 1
        col2 = data[:refnr, colarr[1]]  # Column n       
    return col1, col2
#-------------------------------------------------------------------------------

# Find latest file of a given type
def find_latest_file(flist):
    if not len(flist):
        print(f'ERROR: No file in {flist} found')
        return -1
    outfile = max(flist, key=os.path.getmtime)
    return outfile
#-------------------------------------------------------------------------------
def read_xvg_filedata(fname):
    with open(fname) as fin:
        lines = (
            line.strip()
            for line in fin
            if not line.lstrip().startswith('#')
            and not line.lstrip().startswith('@')
        )

        data = np.loadtxt(lines)

    return data
#-------------------------------------------------------------------------------

# VFT-like model in concentration
# D(c) = D0 * exp[-B / (c0 - c)]
# log(D) = log(D0) - B / (c0 - c)
def log_vft_conc(c, logD0, B, c0):
    return logD0 - (B / (c0 - c))
#-------------------------------------------------------------------------------

def extract_fit_params(filename):
    with open(filename, "r") as f:
        text = f.read()

    data = {}


    block_pattern = re.compile(r"(\w+)\s+log-space fit parameters:\s*(.*?)(?=\n\s*\w+\s+log-space fit parameters:|\Z)",
                               flags=re.IGNORECASE | re.DOTALL)
    
    num = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?"

    # Value and error are separated by comma
    line_pattern = re.compile(rf"([A-Za-z0-9_().]+)\s*=\s*({num})\s*,\s*({num})",
                              flags=re.IGNORECASE)

    for species, block in block_pattern.findall(text):
        data[species] = {}

        for name, value, error in line_pattern.findall(block):
            data[species][name] = {
                "value": float(value),
                "error": float(error)
            }
    return data
#-------------------------------------------------------------------------------

#if __name__
if __name__ == '__main__':
    main()
#-------------------------------------------------------------------------------
