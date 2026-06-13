import MDAnalysis as mda

u = mda.Universe("npt_main.tpr", "traj_npt_main.trr",\
                 refresh_offsets=True)
ts = u.trajectory[-1]
print("n_atoms =", u.atoms.n_atoms)
print("n_frames =", len(u.trajectory))
print("Last time =", ts.time) 
ts = u.trajectory[0]
print("time =", ts.time)
