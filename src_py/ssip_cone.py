# Code to compute the CIP/SSIP and free ion-pairs
# Change the default arguments in aux_ssip.py
# A cone-angle definition is used for SSIP
# Version: May-24-2026
import argparse
import math
import csv
import sys
import numpy as np
import aux_ssip as aip
import MDAnalysis as mda
from MDAnalysis.lib import distances
from MDAnalysis.lib.distances import minimize_vectors


def main():
    args = aip.parse_args()

    u = mda.Universe(args.structure, args.trajectory, refresh_offsets=True)

    print(args.structure,args.trajectory)
    
    fe_atoms = u.select_atoms(args.cation_sel)
    anion_res = u.select_atoms(args.anion_res_sel).residues
    solv_O = u.select_atoms(args.solv_O_sel)

    if len(fe_atoms) == 0:
        raise ValueError("No Fe atoms found with the given --cation-sel")
    
    if len(anion_res) == 0:
        raise ValueError("No TFSI residues found with the given --anion-res-sel")
    
    if len(solv_O) == 0:
        raise ValueError("No water oxygens found with the given --solv-O-sel")
    
    if args.r_tfsi is None:
        r_tfsi, r_tfsi_mean, d_tfsi_mean, n_radius_frames = \
            aip.estimate_tfsi_radius(u,anion_res,tmin=args.tmin,tmax=args.tmax,\
                                     stride=args.stride,mode=args.r_tfsi_mode,\
                                     percentile=args.r_tfsi_percentile
            )
        print(
            f"Estimated TFSI effective radius = {r_tfsi:.3f} Å; "
            f"mean outer-atom radius = {r_tfsi_mean:.3f} Å; "
            f"mean effective diameter = {d_tfsi_mean:.3f} Å "
            f"(mode={args.r_tfsi_mode}, frames_used={n_radius_frames})",
            file=sys.stderr
        )
    else:
        r_tfsi = args.r_tfsi
        _, r_tfsi_mean, d_tfsi_mean, n_radius_frames = \
            aip.estimate_tfsi_radius(u,anion_res,tmin=args.tmin,\
                                     tmax=args.tmax,stride=args.stride,\
                                     mode="percentile",percentile=50.0
            )
        print(
            f"Using user-provided TFSI effective radius = {r_tfsi:.3f} Å; "
            f"mean outer-atom radius = {r_tfsi_mean:.3f} Å; "
            f"mean effective diameter = {d_tfsi_mean:.3f} Å "
            f"(frames_used={n_radius_frames})",
            file=sys.stderr
        )

    labels_history = {}

    ts_out = open(args.out_timeseries, "w", newline="")
    ts_writer = csv.writer(ts_out)
    ts_writer.writerow([
        "frame",
        "time_ps",
        "fe_index",
        "tfsi_resid",
        "label",
        "d_fe_outer_min",
        "d_fe_com",
        "bridged"
    ])

    processed_frames = 0
    last_report_bin = None

    for ts in u.trajectory[::args.stride]:
        tps = float(ts.time)

        if args.tmin is not None and tps < args.tmin:
            continue
        if args.tmax is not None and tps > args.tmax:
            break

        report_bin = int(tps // 1000)
        if report_bin != last_report_bin:
            print(f"Analyzing frame at {tps:.3f} ps")
            last_report_bin = report_bin

        box = u.dimensions

        for fe_atom in fe_atoms.atoms:
            fe_pos = fe_atom.position

            res, dmin = aip.nearest_residue_by_min_outeratom(fe_pos, anion_res, box)
            if res is None:
                continue

            com_pos = aip.residue_com_mic(res, box)
            d_fe_com = np.linalg.norm(aip.mic_vector(fe_pos, com_pos, box))

            if dmin > args.r_assoc:
                label = "free"
                bridged = 0
            else:
                bridge_found = False

                for ow in solv_O.atoms:
                    inside_cone, L, t, theta, alpha = \
                        aip.water_in_tfsi_cone(fe_pos,com_pos,ow.position,\
                                               r_tfsi,box)
                    if inside_cone:
                        bridge_found = True
                        break

                if bridge_found:
                    label = "SSIP"
                    bridged = 1
                else:
                    label = "CIP"
                    bridged = 0

            key = (fe_atom.ix, res.resid)
            labels_history.setdefault(key, []).append(label)

            ts_writer.writerow([ts.frame,f"{tps:.6f}",fe_atom.ix,\
                                res.resid,label,f"{dmin:.3f}",\
                                f"{d_fe_com:.3f}",bridged])

        processed_frames += 1

    ts_out.close()

    cats = ["free", "SSIP", "CIP"]

    raw_counts = aip.tally(labels_history,cats)

    if args.min_residence_frames > 1:
        smoothed_counts = {k: 0 for k in cats}
        for key, seq in labels_history.items():
            sm = aip.smooth_labels(seq, args.min_residence_frames)
            for x in sm:
                if x in smoothed_counts:
                    smoothed_counts[x] += 1
    else:
        smoothed_counts = raw_counts.copy()

    with open(args.out_summary, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["metric"] + cats)
        w.writerow(["tfsi_effective_radius_A", f"{r_tfsi:.6f}"])
        w.writerow(["tfsi_mean_outer_radius_A", f"{r_tfsi_mean:.6f}"])
        w.writerow(["tfsi_mean_effective_diameter_A", f"{d_tfsi_mean:.6f}"])
        w.writerow(["raw_counts"] + [raw_counts[k] for k in cats])

        rf = aip.fracs(raw_counts,cats)
        w.writerow(["raw_fractions"] + [f"{rf[k]:.6f}" for k in cats])

        w.writerow(
            [f"smoothed_counts(n={args.min_residence_frames})"] +
            [smoothed_counts[k] for k in cats]
        )

        sf = aip.fracs(smoothed_counts,cats)
        w.writerow(
            [f"smoothed_fractions(n={args.min_residence_frames})"] +
            [f"{sf[k]:.6f}" for k in cats]
        )

    print(f"Processed {processed_frames} frame(s).")
    print(f"Wrote {args.out_timeseries} and {args.out_summary}")


if __name__ == "__main__":
    main()
