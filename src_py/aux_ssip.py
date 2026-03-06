"""
Auxiliary files for ssip.py
"""
import argparse, math, csv, sys
import numpy as np
import MDAnalysis as mda
from MDAnalysis.lib import distances

DEG = math.pi/180.0
#-----------------------------------------------------------------------------------------
def unit(v):
    n = np.linalg.norm(v)
    if n == 0: return v
    return v / n
#-----------------------------------------------------------------------------------------

def angle(a, b, c):
    """Return angle at vertex b (∠a-b-c) in degrees."""
    v1 = a - b; v2 = c - b
    u1, u2 = unit(v1), unit(v2)
    x = np.clip(np.dot(u1, u2), -1.0, 1.0)
    return math.degrees(math.acos(x))
#-----------------------------------------------------------------------------------------

def min_dist_atom_to_group(atom_pos, group_pos, box):
    """Minimum-image distance from one position to ANY in group."""
    d = distances.distance_array(atom_pos[np.newaxis,:], group_pos, box=box, backend='OpenMP')
    return float(d.min())
#-----------------------------------------------------------------------------------------

def nearest_residue_COM(atom_pos, res_COMs, box):
    d = distances.distance_array(atom_pos[np.newaxis,:], res_COMs, box=box, backend='OpenMP')
    idx = int(np.argmin(d[0]))
    return idx, float(d[0, idx])
#-----------------------------------------------------------------------------------------

def project_geometry(Fe, AnionCOM, S, eps=1e-8):
    """
    Return (L, t, d_line, phi_deg) with:
      v = AnionCOM - Fe ; L = |v| ; u = v/|v|
      w = S - Fe ; t = w·u ; d_line^2 = |w|^2 - t^2 ; phi = asin(d_line/|w|)
    All in Cartesian coords (assumes unwrapped or single image continuity).
    """
    v = AnionCOM - Fe
    L = np.linalg.norm(v)
    if L < eps:
        return 0.0, 0.0, 0.0, 0.0
    
    u = v / L
    w = S - Fe
    t = float(np.dot(w, u))
    perp2 = float(np.dot(w, w) - t*t)
    if perp2 < 0.0:  # numerical safety
        perp2 = 0.0
    d_line = math.sqrt(perp2)
    wnorm = np.linalg.norm(w)
    phi = 0.0 if wnorm < eps else math.degrees(math.asin(max(0.0, min(1.0, d_line/max(wnorm, eps)))))
    return L, t, d_line, phi
#-----------------------------------------------------------------------------------------

def smooth_labels(seq, min_len):
    """Hysteresis smoothing: convert short runs (<min_len) to neighboring label."""
    if min_len <= 1 or not seq:
        return seq[:]
    out = []
    i = 0
    n = len(seq)
    while i < n:
        j = i
        while j < n and seq[j] == seq[i]:
            j += 1
        run_len = j - i
        if run_len < min_len:
            fill = out[-1] if out else (seq[j] if j < n else seq[i])
            out.extend([fill]*run_len)
        else:
            out.extend([seq[i]]*run_len)
        i = j
    return out
#-----------------------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    p.add_argument('-s','--structure', required=True, help='Topology (e.g., .tpr/.gro/.pdb)')
    p.add_argument('-f','--trajectory', required=True, help='Trajectory (e.g., .trr/.xtc/.dcd)')
    # selections
    p.add_argument('--cation-sel', default='name Fe or name FE or element Fe', help='Fe cation selection')
    p.add_argument('--anion-res-sel', default='resname TFSI or resname TFS or resname T2N', help='Residues comprising TFSI anions')
    p.add_argument('--anion-O-sel', default='(resname TFSI or resname TFS or resname T2N) and name O*', help='O atoms of TFSI')
    p.add_argument('--solv-res-sel', default='resname SOL or resname WAT', help='Solvent (water) residues')
    p.add_argument('--solv-O-sel', default='(resname SOL and name OW) or (resname WAT and name O)', help='Water oxygen atoms')
    # Iteration / time window
    p.add_argument('--stride', type=int, default=10, help='Process every Nth frame')
    p.add_argument('--tmin', type=float, default=90000, help='Start time (ps), inclusive')
    p.add_argument('--tmax', type=float, default=92000, help='End time (ps), inclusive')
    # geometry params
    p.add_argument('--tube-radius', type=float, default=1.8, help='Tube radius Å (size-aware default for Ow)')
    p.add_argument('--phi-max', type=float, default=30.0, help='Max cone half-angle in degrees')
    # chemistry params
    p.add_argument('--r-feO-contact', type=float, default=2.6,
                   help='Fe–O_TFSI contact cutoff (Å), typical ~2.5–2.7; REQUIRED for robust CIP')
    p.add_argument('--r-feOw', type=float, default=2.2,
                   help='Fe–O_w first-shell cutoff (Å), typical ~2.3–2.6; REQUIRED for strict SSIP')
    p.add_argument('--hbond-OO', type=float, default=3.5, help='Ow···O_TFSI O–O cutoff for H-bond (Å)')
    p.add_argument('--hbond-angle', type=float, default=150.0, help='∠H–Ow···O_TFSI (deg) minimum')
    # smoothing
    p.add_argument('--min-residence-frames', type=int, default=5, help='Time smoothing (frames)')
    # output
    p.add_argument('--out-timeseries', default='pairs_timeseries.csv')
    p.add_argument('--out-summary', default='pairs_summary.csv')
    return p.parse_args()
#-----------------------------------------------------------------------------------------
def tally(seqs,cats):
    counts = {k:0 for k in cats}
    for seq in seqs.values():
        for L in seq:
            counts[L] += 1
    return counts
#-----------------------------------------------------------------------------------------
def to_frac(counts,cats):
    tot = sum(counts.values())
    return {k:(counts[k]/tot if tot else 0.0) for k in cats}
