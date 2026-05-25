import argparse
import math
import csv
import sys
import numpy as np
import MDAnalysis as mda
from MDAnalysis.lib import distances

def unit(v):
    n = np.linalg.norm(v)
    return v if n == 0 else v / n

def angle(a, b, c):
    """Return angle at vertex b (∠a-b-c) in degrees."""
    v1, v2 = a - b, c - b
    u1, u2 = unit(v1), unit(v2)
    x = np.clip(np.dot(u1, u2), -1.0, 1.0)
    return math.degrees(math.acos(x))

def project_geometry(Fe, AnionCOM, S, eps=1e-8):
    """
    Geometry for already-unwrapped coordinates.
    Inputs must be shape (3,).
    Returns:
      L      = |Fe -> AnionCOM|
      t      = projection of (Fe -> S) along Fe -> AnionCOM
      d_line = perpendicular distance of S from Fe-Anion axis
      phi    = angle between (Fe -> S) and (Fe -> AnionCOM), degrees
    """
    Fe = np.asarray(Fe, dtype=float).reshape(3,)
    AnionCOM = np.asarray(AnionCOM, dtype=float).reshape(3,)
    S = np.asarray(S, dtype=float).reshape(3,)

    v = AnionCOM - Fe
    L = np.linalg.norm(v)
    if L < eps:
        return 0.0, 0.0, 0.0, 0.0

    u = v / L
    w = S - Fe

    t = float(np.dot(w, u))

    perp2 = float(np.dot(w, w) - t * t)
    if perp2 < 0.0:
        perp2 = 0.0
    d_line = math.sqrt(perp2)

    wnorm = np.linalg.norm(w)
    phi = 0.0 if wnorm < eps else math.degrees(
        math.asin(np.clip(d_line / max(wnorm, eps), 0.0, 1.0))
    )

    return L, t, d_line, phi

def min_dist_atom_to_group(atom_pos, group_pos, box):
    d = distances.distance_array(atom_pos[np.newaxis, :], group_pos, box=box, backend='OpenMP')
    return float(d.min())

def nearest_residue_COM(atom_pos, res_COMs, box):
    d = distances.distance_array(atom_pos[np.newaxis, :], res_COMs, box=box, backend='OpenMP')
    idx = int(np.argmin(d[0]))
    return idx, float(d[0, idx])

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

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('-s', '--structure', required=True)
    p.add_argument('-f', '--trajectory', required=True)

    # selections
    p.add_argument('--cation-sel', default='name Fe or name FE or element Fe', help='Fe cation selection')
    p.add_argument('--anion-res-sel', default='resname TFSI or resname TFS or resname T2N', help='Residues comprising TFSI anions')
    p.add_argument('--anion-O-sel', default='(resname TFSI or resname TFS or resname T2N) and type OS', help='O atoms of TFSI')
    p.add_argument('--solv-res-sel', default='resname SOL or resname WAT', help='Solvent (water) residues')
    p.add_argument('--solv-O-sel', default='(resname SOL and name OW) or (resname WAT and name O)', help='Water oxygen atoms')

    # RDF-derived shell cutoffs
    p.add_argument('--r-feOw', type=float, default=2.4,
                   help='Fe-Ow first-shell cutoff in Å')
    p.add_argument('--r-assoc', type=float, default=5.3,
                   help='Fe-OTFSI association-shell cutoff in Å')

    # bridge geometry
    p.add_argument('--tube-radius', type=float, default=1.8)
    p.add_argument('--phi-max', type=float, default=30.0)

    # H-bond-like Ow...OS bridge criterion
    p.add_argument('--hbond-OO', type=float, default=3.5,
                   help='Ow...OS cutoff in Å')
    p.add_argument('--hbond-angle', type=float, default=150.0,
                   help='H-Ow...OS minimum angle in degrees')

    # Iteration / time window
    p.add_argument('--stride', type=int, default=100, help='Process every Nth frame')
    p.add_argument('--tmin', type=float, default=80000, help='Start time (ps), inclusive')
    p.add_argument('--tmax', type=float, default=90000, help='End time (ps), inclusive')
    p.add_argument('--min-residence-frames', type=int, default=10, help='Time smoothing (frames)')

    p.add_argument('--out-timeseries', default='newpairs_timeseries.csv')
    p.add_argument('--out-summary', default='newpairs_summary.csv')
    return p.parse_args()

def main():
    args = parse_args()
    print(args.structure, args.trajectory)

    u = mda.Universe(args.structure, args.trajectory, refresh_offsets=True)

    Fe = u.select_atoms(args.cation_sel)
    anion_res = u.select_atoms(args.anion_res_sel).residues
    anion_O_all = u.select_atoms(args.anion_O_sel)
    solv_res = u.select_atoms(args.solv_res_sel).residues
    solv_O = u.select_atoms(args.solv_O_sel)

    if len(Fe) == 0 or len(anion_res) == 0 or len(solv_O) == 0:
        raise ValueError("One or more atom selections are empty.")

    # Precompute solvent Hs per water residue (used for H-bond angle)
    H_by_solvent_resid = {}
    for res in solv_res:
        H_by_solvent_resid[res.resid] = res.atoms.select_atoms('element H or name H*')

    labels_history = {}

    ts_out = open(args.out_timeseries, 'w', newline='')
    ts_writer = csv.writer(ts_out)
    ts_writer.writerow([
        'frame', 'time_ps', 'fe_index', 'tfsi_resid', 'label',
        'd_feO_tfsi_min', 'd_fe_anionCOM', 'bridged'
    ])

    nframes = len(u.trajectory)
    processed_frames = 0

    # Iterate frames
    for ts in u.trajectory[::args.stride]:
        tps = float(ts.time)

        if args.tmin is not None and tps < args.tmin:
            continue
        if args.tmax is not None and tps > args.tmax:
            break


        if args.tmin is not None and (tps-args.tmin)%100 == 0:
            print(f"Analyzing frame at {tps} ps")

        box = u.dimensions
        anion_COMs = np.array([res.atoms.center_of_mass() for res in anion_res], dtype=float)

        
        for fe_atom in Fe.atoms:
            fe_pos = fe_atom.position

            idx, dCOM = nearest_residue_COM(fe_pos, anion_COMs, box)
            res = anion_res[idx]
            resCOM = anion_COMs[idx]
            resO = res.atoms.select_atoms('type OS')

            dmin = min_dist_atom_to_group(fe_pos, resO.positions, box) if len(resO) else float('inf')

            # Outside association shell
            if dmin > args.r_assoc:
                label = 'free'
                bridged = 0

            else:
                bridge_found = False

                for ow in solv_O.atoms:
                    ow_pos = ow.position

                    L, t, d_line, phi = project_geometry(fe_pos, resCOM, ow_pos)

                    if L <= 1e-8:
                        continue
                    if not (t > 1e-6 and t < L - 1e-6):
                        continue
                    if d_line > args.tube_radius:
                        continue
                    if phi > args.phi_max:
                        continue

                    d_FeOw = distances.distance_array(
                        fe_pos[np.newaxis, :], ow_pos[np.newaxis, :], box=box, backend='OpenMP'
                    )[0, 0]
                    if d_FeOw > args.r_feOw:
                        continue

                    hbond_ok = False
                    dOO_all = distances.distance_array(
                        ow_pos[np.newaxis, :], resO.positions, box=box, backend='OpenMP'
                    )[0]
                    close_idxs = np.where(dOO_all <= args.hbond_OO)[0]

                    if close_idxs.size > 0:
                        solvent_Hs = H_by_solvent_resid.get(ow.resid, u.atoms[:0])
                        for oi in close_idxs:
                            os_pos = resO.positions[oi]
                            for H in solvent_Hs:
                                ang = angle(H.position, ow_pos, os_pos)
                                if ang >= args.hbond_angle:
                                    hbond_ok = True
                                    break
                            if hbond_ok:
                                break

                    if hbond_ok:
                        bridge_found = True
                        break

                if bridge_found:
                    label = 'SSIP'
                    bridged = 1
                else:
                    label = 'CIP'
                    bridged = 0

            key = (fe_atom.ix, res.resid)
            labels_history.setdefault(key, []).append(label)

            ts_writer.writerow([
                ts.frame, f"{tps:.6f}", fe_atom.ix, res.resid, label,
                f"{dmin:.3f}", f"{dCOM:.3f}", bridged
            ])

        processed_frames += 1

    ts_out.close()

    cats = ['free', 'SSIP', 'CIP']

    def tally(hist):
        counts = {k: 0 for k in cats}
        for seq in hist.values():
            for x in seq:
                if x in counts:
                    counts[x] += 1
        return counts

    def fracs(counts):
        total = sum(counts.values())
        return {k: (counts[k] / total if total else 0.0) for k in cats}

    raw_counts = tally(labels_history)

    smoothed_counts = {k: 0 for k in cats}
    if args.min_residence_frames > 1:
        for key, seq in labels_history.items():
            sm = smooth_labels(seq, args.min_residence_frames)
            for x in sm:
                if x in smoothed_counts:
                    smoothed_counts[x] += 1
    else:
        smoothed_counts = raw_counts.copy()

    with open(args.out_summary, 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['metric'] + cats)
        w.writerow(['raw_counts'] + [raw_counts[k] for k in cats])
        rf = fracs(raw_counts)
        w.writerow(['raw_fractions'] + [f"{rf[k]:.6f}" for k in cats])
        w.writerow([f'smoothed_counts(n={args.min_residence_frames})'] + [smoothed_counts[k] for k in cats])
        sf = fracs(smoothed_counts)
        w.writerow([f'smoothed_fractions(n={args.min_residence_frames})'] + [f"{sf[k]:.6f}" for k in cats])

    print(f"Processed {processed_frames} frame(s).")
    print(f"Wrote {args.out_timeseries} and {args.out_summary}")

if __name__ == '__main__':
    main()
