"""
Auxiliary files for ssip.py
"""
import argparse, math, csv, sys
import numpy as np
import MDAnalysis as mda
from MDAnalysis.lib import distances
from MDAnalysis.lib.distances import minimize_vectors

DEG = math.pi/180.0
#-----------------------------------------------------------------------------------------
def unit(v):
    n = np.linalg.norm(v)
    return v if n == 0 else v / n
#-----------------------------------------------------------------------------------------
def mic_vector(r1, r2, box):
    v = np.asarray(r2, dtype=float).reshape(1, 3) - np.asarray(r1, dtype=float).reshape(1, 3)
    return minimize_vectors(v, box)[0]
#-----------------------------------------------------------------------------------------
def unwrap_group_positions(group, box):
    pos = group.positions.copy()
    ref = pos[0].copy()
    rel = pos - ref
    rel = minimize_vectors(rel, box)
    return ref + rel
#-----------------------------------------------------------------------------------------
def residue_com_mic(res, box):
    pos_unwrapped = unwrap_group_positions(res.atoms, box)
    masses = res.atoms.masses
    return np.sum(pos_unwrapped * masses[:, None], axis=0) / np.sum(masses)
#-----------------------------------------------------------------------------------------
def min_dist_atom_to_group(atom_pos, group_pos, box):
    d = distances.distance_array(
        atom_pos[np.newaxis, :],
        group_pos,
        box=box,
        backend="OpenMP"
    )
    return float(d.min())
#-----------------------------------------------------------------------------------------
def nearest_residue_by_min_outeratom(fe_pos, anion_res, box):
    """
    Choose the TFSI residue whose relevant outer atoms (OS/FC)
    have the minimum Fe distance.
    """
    best_res = None
    best_dmin = float("inf")

    for res in anion_res:
        outer_atoms = res.atoms.select_atoms("type OS or type FC")
        if outer_atoms.n_atoms == 0:
            continue

        dmin = min_dist_atom_to_group(fe_pos, outer_atoms.positions, box)
        if dmin < best_dmin:
            best_dmin = dmin
            best_res = res

    return best_res, best_dmin
#-----------------------------------------------------------------------------------------
def water_in_tfsi_cone(fe_pos,com_pos,ow_pos,r_tfsi,box,eps=1e-8):
    """
    Use a cone to separate the Fe-TFSI and water atom
    A water is counted as separating if its direction from Fe lies
    inside the cone subtended by the sphere and it lies between Fe
    & the near surface of the sphere 
    """
    fe_pos = np.asarray(fe_pos, dtype=float).reshape(3,)
    com_pos = np.asarray(com_pos, dtype=float).reshape(3,)
    ow_pos = np.asarray(ow_pos, dtype=float).reshape(3,)

    r_fc = mic_vector(fe_pos, com_pos, box)
    L = np.linalg.norm(r_fc)
    if L < eps:
        return False, 0.0, 0.0, 0.0, 0.0

    ratio = min(max(r_tfsi / L, 0.0), 1.0)
    alpha = math.degrees(math.asin(ratio))

    u_fc = r_fc / L
    r_fw = mic_vector(fe_pos, ow_pos, box)
    rw = np.linalg.norm(r_fw)
    if rw < eps:
        return False, L, 0.0, 0.0, alpha

    t = float(np.dot(r_fw, u_fc))

    cos_theta = np.clip(np.dot(r_fw, r_fc) / (rw * L), -1.0, 1.0)
    theta = math.degrees(math.acos(cos_theta))

    inside_cone = (theta < alpha) and (0.0 < t < (L - r_tfsi))

    return inside_cone, L, t, theta, alpha
#-----------------------------------------------------------------------------------------
def estimate_tfsi_radius(u,anion_res,tmin=None,tmax=None,stride=1,\
                         mode="percentile",percentile=95.0):
    """
    Estimating TFSI radius is an approximation. Here, the effective
    TFSI radius from distances of outer atoms (OS/FC) from the residue
    COM over sampled frames since Fe attaches to both OS and FC atoms
    """

    nframes = len(u.trajectory)
    if nframes == 0:
        raise ValueError("Trajectory has no frames; cannot estimate TFSI radius.")

    radii = []
    n_used_frames = 0
    current_frame = u.trajectory.frame

    for ts in u.trajectory[::stride]:
        tps = float(ts.time)

        if tmin is not None and tps < tmin:
            continue
        if tmax is not None and tps > tmax:
            break

        box = u.dimensions
        n_used_frames += 1

        for res in anion_res:
            outer_atoms = res.atoms.select_atoms("type OS or type FC")
            if outer_atoms.n_atoms == 0:
                continue

            com = residue_com_mic(res, box)

            # MIC-consistent displacement of outer atoms from COM
            rel = outer_atoms.positions - com
            rel = minimize_vectors(rel, box)
            d = np.linalg.norm(rel, axis=1)
            radii.extend(d.tolist())

    _ = u.trajectory[current_frame]

    if n_used_frames == 0:
        raise ValueError("No frames in the requested tmin/tmax/stride window for TFSI radius estimation.")

    if len(radii) == 0:
        raise ValueError("Could not estimate TFSI radius: no OS/FC atoms found.")

    radii = np.asarray(radii, dtype=float)

    r_mean = float(np.mean(radii))
    d_mean = 2.0 * r_mean

    if mode == "max":
        r_eff = float(np.max(radii))
    else:
        r_eff = float(np.percentile(radii, percentile))

    return r_eff, r_mean, d_mean, n_used_frames
#-----------------------------------------------------------------------------------------
def smooth_labels(seq, min_len):
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
            out.extend([fill] * run_len)
        else:
            out.extend([seq[i]] * run_len)

        i = j

    return out
#-----------------------------------------------------------------------------------------
def parse_args():
    p = argparse.ArgumentParser()

    # Trajectory parsers
    p.add_argument("-s", "--structure", required=True)
    p.add_argument("-f", "--trajectory", required=True)

    # Selections
    p.add_argument('--cation-sel', default='name Fe or name FE or element Fe', help='Fe cation selection')
    p.add_argument('--anion-res-sel', default='resname TFSI or resname TFS or resname T2N', help='Residues comprising TFSI anions')
    p.add_argument('--anion-O-sel', default='(resname TFSI or resname TFS or resname T2N) and type OS', help='O atoms of TFSI')
    p.add_argument('--solv-res-sel', default='resname SOL or resname WAT', help='Solvent (water) residues')
    p.add_argument('--solv-O-sel', default='(resname SOL and name OW) or (resname WAT and name O)', help='Water oxygen atoms')
    
    # Geometry constraints
    p.add_argument("--r-assoc",type=float,default=5.3,help="Fe-TFSI association-shell cutoff in Å, based on outer atoms (OS/FC)"
    )

    p.add_argument("--r-tfsi",type=float,default=None,help="Effective TFSI sphere radius in Å. If omitted, estimate from trajectory."
    )

    p.add_argument("--r-tfsi-mode",choices=["percentile", "max"],\
                   default="percentile",help="How to estimate TFSI radius from OS/FC atoms if --r-tfsi is omitted"
    )

    p.add_argument("--r-tfsi-percentile",type=float,default=95.0,
                   help="Percentile for TFSI radius estimation when --r-tfsi-mode=percentile"
    )

    # I/O arguments
    p.add_argument("--stride", type=int, default=500)
    p.add_argument("--tmin", type=float, default=80000)
    p.add_argument("--tmax", type=float, default=90000)
    p.add_argument("--min-residence-frames", type=int, default=10)

    p.add_argument("--out-timeseries", default="conepairs_timeseries.csv")
    p.add_argument("--out-summary", default="conepairs_summary.csv")

    return p.parse_args()
#-----------------------------------------------------------------------------------------
def angle(a, b, c):
    """Return angle at vertex b (∠a-b-c) in degrees."""
    v1 = a - b; v2 = c - b
    u1, u2 = unit(v1), unit(v2)
    x = np.clip(np.dot(u1, u2), -1.0, 1.0)
    return math.degrees(math.acos(x))
#-----------------------------------------------------------------------------------------
def tally(hist,cats):
    counts = {k: 0 for k in cats}
    for seq in hist.values():
        for x in seq:
            if x in counts:
                counts[x] += 1
    return counts
#-----------------------------------------------------------------------------------------
def fracs(counts,cats):
    total = sum(counts.values())
    return {k: (counts[k] / total if total else 0.0) for k in cats}
#-----------------------------------------------------------------------------------------
def nearest_residue_COM(atom_pos, res_COMs, box):
    d = distances.distance_array(atom_pos[np.newaxis,:], res_COMs, box=box, backend='OpenMP')
    idx = int(np.argmin(d[0]))
    return idx, float(d[0, idx])
#-----------------------------------------------------------------------------------------
def project_geometry(Fe, AnionCOM, OW_pos, box, eps=1e-8):
    """
    Return (L, t, d_line, phi_deg) with:
      v = AnionCOM - Fe ; L = |v| ; u = v/|v|
      w = OW_pos - Fe ; t = w·u ; d_line^2 = |w|^2 - t^2 ; phi = asin(d_line/|w|)
    All in Cartesian coords (assumes unwrapped or single image continuity).
    """

    v = mic_vector(Fe,AnionCOM,box)
    L = np.linalg.norm(v)

    if L < eps:
        return 0.0, 0.0, 0.0, 0.0
    u = v / L

    w = mic_vector(Fe,OW_pos,box)
    
    t = float(np.dot(w, u))
    perp2 = float(np.dot(w, w) - t*t)

    if perp2 < 0.0:  # numerical safety
        perp2 = 0.0
    d_line = math.sqrt(perp2)

    wnorm = np.linalg.norm(w)
    phi = 0.0 if wnorm < eps else np.degrees(
        np.arcsin(np.clip(d_line / max(wnorm, eps), 0.0, 1.0))
    )
    
    return L, t, d_line, phi
#-----------------------------------------------------------------------------------------
