# check that the length of dtrajs in order matches with the length of the trajectories in order
 
import numpy as np
# discrete trajectories first
def load_dtrajs_from_npz(filename):
    dtrajs = [y[1] for y in np.load(filename).items()]
    return dtrajs

dtrajs_abl = load_dtrajs_from_npz('../results/abl_11400_2000_dtrajs.npz')
dtrajs_src = load_dtrajs_from_npz('../results/src_11401_2000_dtrajs.npz')

dtraj_lengths_abl = [len(t) for t in dtrajs_abl]
dtraj_lengths_src = [len(t) for t in dtrajs_src]

# now actual trajectories
from glob import glob
traj_files_abl = glob('../abl_snapshot/*.h5')
traj_files_src = glob('../src_snapshot/*.h5')

import mdtraj as md
traj_lengths_abl = [len(md.load(t)) for t in traj_files_abl]
traj_lengths_src = [len(md.load(t)) for t in traj_files_src]

# print comparison
print('Abl: (dtraj_lengths,traj_lengths)')
print(zip(dtraj_lengths_abl,traj_lengths_abl))

print('Src: (dtraj_lengths,traj_lengths)')
print(zip(dtraj_lengths_src,traj_lengths_src))
