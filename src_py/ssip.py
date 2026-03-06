#!/usr/bin/env python3
"""
Classify Fe–TFSI ion pairs in water as CIP vs SSIP (size-aware + chemistry-aware).

Definitions (per frame, nearest TFSI per Fe):
- CIP:    min( Fe–O_TFSI ) <= r_FeO_contact  (contact by Fe–O_TFSI first-min from g(r))
- SSIP:   NOT CIP AND exists water oxygen (Ow) that is:
          (i)  geometrically BETWEEN Fe and the TFSI COM:
               0 < t = (Ow-Fe)·u < L  AND
               d_line = distance from Ow to Fe–TFSI axis <= R_tube  AND
               cone half-angle φ <= φ_max
          (ii) chemistry-aware bridge:
               Fe–Ow <= r_FeOw  AND  H-bond Ow···O_TFSI (O–O <= 3.5 Å and ∠H–Ow···O_TFSI >= 150°)
- SSIP_loose: passes geometry but not chemistry (optional label)

Outputs:
  - timeseries CSV: frame, time_ps, fe_ix, tfsi_resid, label, d_fe_otfsi_min, d_fe_anionCOM, bridged(0/1)
  - summary CSV: fractions of CIP/SSIP/SSIP_loose/unassigned (raw and smoothed)
"""

import argparse, math, csv, sys
import numpy as np
import MDAnalysis as mda
from MDAnalysis.lib import distances
import aux_ssip as ssp
from MDAnalysis.transformations import unwrap

DEG = math.pi/180.0
#-----------------------------------------------------------------------------------------
def main():
    args = ssp.parse_args()
    print(f'Structure-file: {args.structure}; Trajectory-file: {args.trajectory}')
    u = mda.Universe(args.structure, args.trajectory,refresh_offsets=True)

#    Try to unwrap residues (so geometry is continuous under PBC)
    try:
        u.trajectory.add_transformations(unwrap(u.atoms))
    except Exception:
        print("Note: could not enable on-the-fly unwrap; ensure your box is orthorhombic or pre-unwrapped.", file=sys.stderr)
 
    # Fe = u.select_atoms(args.cation-sel)  # will fail; hyphen not allowed in attr access
    # Fix: use getattr style
    Fe = u.select_atoms(args.cation_sel)
    anion_res = u.select_atoms(args.anion_res_sel).residues
    anion_O = u.select_atoms(args.anion_O_sel)
    solv_res = u.select_atoms(args.solv_res_sel).residues
    solv_O = u.select_atoms(args.solv_O_sel)

    if len(Fe) == 0 or len(anion_res) == 0 or len(solv_O) == 0:
        raise ValueError("Empty selection: check your --*_sel strings.")

    if args.r_feO_contact is None or args.r_feOw is None:
        print("WARNING: --r-feO-contact and --r-feOw not provided; "
              "set them from your RDF first minima (Fe–O_TFSI and Fe–Ow). "
              "Using heuristics may misclassify.", file=sys.stderr)
    r_FeO_contact = args.r_feO_contact if args.r_feO_contact is not None else 2.6
    r_FeOw = args.r_feOw if args.r_feOw is not None else 2.4

    # Precompute atoms per residue for TFSI O atoms
    tfsi_O_by_resid = {}
    for res in anion_res:
        tfsi_O_by_resid[res.resid] = res.atoms.select_atoms('name O*') if res.atoms.n_atoms > 0 else u.atoms[:0]

    # Precompute solvent Hs per water residue (used for H-bond angle)
    H_by_solvent_resid = {}
    for res in solv_res:
        H_by_solvent_resid[res.resid] = res.atoms.select_atoms('element H or name H*')

    # Output writers
    ts_out = open(args.out_timeseries, 'w', newline='')
    ts_writer = csv.writer(ts_out)
    ts_writer.writerow(['frame','time_ps','fe_index','tfsi_resid','label','d_feO_tfsi_min','d_fe_anionCOM','bridged'])

    # Store labels over time for smoothing
    labels_history = {}  # key: (fe_ix, tfsi_resid) -> [labels]

    # Start analysis here
    processed_frames = 0
    
    # Iterate frames
    for ts in u.trajectory[::args.stride]:
        tps = float(ts.time)

        if tps > args.tmax:
            print(f'Exceeded  {args.tmax}')
            break
        elif (args.tmin is not None and tps < args.tmin):
            continue

        print(f'Starting analysis at {tps}')
        box = u.dimensions  # [lx, ly, lz, alpha, beta, gamma]
        # COM of each TFSI residue
        anion_COMs = np.array([res.atoms.center_of_mass() for res in anion_res], dtype=float)

        # For each Fe, find nearest TFSI (by COM), then classify CIP/SSIP
        for fe_atom in Fe.atoms:
            fe_pos = fe_atom.position
            idx, dCOM = ssp.nearest_residue_COM(fe_pos, anion_COMs, box)
            res = anion_res[idx]
            resCOM = anion_COMs[idx]
            resO = tfsi_O_by_resid.get(res.resid, res.atoms.select_atoms('name O*'))

            # CIP: any Fe–O_TFSI within contact cutoff?
            dmin = ssp.min_dist_atom_to_group(fe_pos, resO.positions, box) if len(resO) else float('inf')
            if dmin <= r_FeO_contact:
                label = 'CIP'
                bridged = 0
            else:
                # Geometric "between": scan candidate water oxygens
                label = 'unassigned'
                bridged = 0
                # Optional quick filter: Ow close to either Fe or anion COM
                # (limits work for big boxes)
                # Keep all if you prefer maximum sensitivity.
                for ow in solv_O.atoms:
                    ow_pos = ow.position
                    L, t, d_line, phi = ssp.project_geometry(fe_pos, resCOM, ow_pos)
                    if L <= 1e-6:   # degenerate
                        continue
                    # strictly interior and inside tube + cone
                    if not (t > 1e-6 and t < L - 1e-6):
                        continue
                    if d_line > args.tube_radius:
                        continue
                    if phi > args.phi_max:
                        continue

                    # Chemistry checks: Fe–Ow and H-bond to some O_TFSI
                    d_FeOw = distances.distance_array(fe_pos[np.newaxis,:], ow_pos[np.newaxis,:], box=box)[0,0]
                    if d_FeOw > r_FeOw:
                        continue
                    # H-bond (Ow···O_TFSI)
                    hbond_ok = False
                    if len(resO) > 0:
                        # O–O distance screen
                        dOO = distances.distance_array(ow_pos[np.newaxis,:], resO.positions, box=box)
                        close_mask = (dOO[0] <= args.hbond_OO)
                        if np.any(close_mask):
                            # angle H–Ow···O_TFSI
                            solvent_Hs = H_by_solvent_resid.get(ow.resid, u.atoms[:0])
                            if len(solvent_Hs) > 0:
                                close_Oidxs = np.where(close_mask[0])[0]
                                for oi in close_Oidxs:
                                    Ot_pos = resO.positions[oi]
                                    for H in solvent_Hs:
                                        ang = ssp.angle(H.position, ow_pos, Ot_pos)
                                        if ang >= args.hbond_angle:
                                            hbond_ok = True
                                            break
                                    if hbond_ok:
                                        break
                    if hbond_ok:
                        label = 'SSIP'
                        bridged = 1
                        break
                    else:
                        # purely geometric bridge; keep as a fallback label
                        label = 'SSIP_loose'
                        bridged = 1
                        # do NOT break; keep searching for a chemistry-validated bridge
                # end loop over Ow
            # record
            key = (fe_atom.ix, res.resid)
            labels_history.setdefault(key, []).append(label)
            ts_writer.writerow([ts.frame,f"{tps:.6f}" , fe_atom.ix, res.resid,\
                                label, f"{dmin:.3f}", f"{dCOM:.3f}",\
                                bridged])


        processed_frames += 1
        
    ts_out.close()

    if processed_frames == 0:
        raise RuntimeError("WARNING: No frames were processed inside the requested time window; "
                           "check --tmin/--tmax against your trajectory times (ps).", file=sys.stderr)
    
    # Summaries (raw and smoothed)
    cats = ['CIP','SSIP','SSIP_loose','unassigned']
    raw_counts = ssp.tally(labels_history,cats)
    smoothed_counts = {k:0 for k in cats}

    if args.min_residence_frames > 1:
        for key, seq in labels_history.items():
            sm = ssp.smooth_labels(seq, args.min_residence_frames)
            for L in sm:
                smoothed_counts[L] += 1
    else:
        smoothed_counts = raw_counts.copy()

    with open(args.out_summary, 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['metric'] + cats)
        w.writerow(['raw_counts'] + [raw_counts[k] for k in cats])
        rf = ssp.to_frac(raw_counts,cats)
        w.writerow(['raw_fractions'] + [f"{rf[k]:.6f}" for k in cats])
        w.writerow([f'smoothed_counts(n={args.min_residence_frames})'] + [smoothed_counts[k] for k in cats])
        sf = ssp.to_frac(smoothed_counts,cats)
        w.writerow([f'smoothed_fractions(n={args.min_residence_frames})'] + [f"{sf[k]:.6f}" for k in cats])

    # Friendly footer
    twin = f"[{args.tmin if args.tmin is not None else '-inf'}, {args.tmax if args.tmax is not None else '+inf'}] ps"
    print(f"Processed {processed_frames} frame(s) within time window {twin}.")
    print(f"Wrote: {args.out_timeseries} and {args.out_summary}")
    print(f"Wrote: {args.out_timeseries} and {args.out_summary}")
#-----------------------------------------------------------------------------------------
if __name__ == "__main__":
    main()
