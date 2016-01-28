# samples from non-coarse-grained MSM

import numpy as np
import numpy.random as npr
import mdtraj as md

# load files
def load_dtrajs_from_npz(filename):
    dtrajs = [y[1] for y in np.load(filename).items()]
    return dtrajs

# these will have the same lengths as the source trajectory files
raw_dtrajs_abl = load_dtrajs_from_npz('../results/abl_11400_2000_dtrajs.npz')
raw_dtrajs_src = load_dtrajs_from_npz('../results/src_11401_2000_dtrajs.npz')

def trim_dtrajs(dtrajs,length_to_discard=10,min_length=100):
   ''' Ignore trajectories under a specified length, and discard
   the initial frames of the remainder.

   Parameters
   ----------
   dtrajs : list of arrays

   length_to_discard : int, optional
   min_length : int, optional

   Returns
   -------
   trimmed : list of arrays
   '''
   trimmed = [traj[length_to_discard:] for traj in dtrajs if len(traj)>length_to_discard]
   return trimmed

# for MSM estimation, ignore initial bit
dtrajs_abl = trim_dtrajs(raw_dtrajs_abl)
dtrajs_src = trim_dtrajs(raw_dtrajs_src)

# estimate MSMs
import msmbuilder.msm

msm_abl = msmbuilder.msm.MarkovStateModel(lag_time=40,ergodic_cutoff=1)
msm_abl.fit(dtrajs_abl)
f = open('abl_msm_report.txt','w')
f.writelines(msm_abl.summarize())
f.close()

msm_src = msmbuilder.msm.MarkovStateModel(lag_time=40,ergodic_cutoff=1)
msm_src.fit(dtrajs_src)
f = open('src_msm_report.txt','w')
f.writelines(msm_src.summarize())
f.close()

# draw samples
n_samples_per_state=1
abl_samples = msm_abl.draw_samples(raw_dtrajs_abl,n_samples_per_state)
src_samples = msm_src.draw_samples(raw_dtrajs_abl,n_samples_per_state)

# # coarse-grain abl?
# from msmbuilder import lumping
# n_macrostates = 16
# pcca_abl = lumping.PCCAPlus.from_msm(msm_abl, n_macrostates=n_macrostates)
# abl_samples = pcca.draw_samples(dtrajs_abl,n_samples_per_state)
#
# # coarse-grain src?
# n_macrostates = 24
# pcca = lumping.PCCAPlus.from_msm(msm_src, n_macrostates=n_macrostates)
# src_samples = pcca.draw_samples(dtrajs_abl,n_samples_per_state)

# fetch corresponding configurations
from glob import glob
traj_files_abl = glob('../abl_snapshot/*.h5')
traj_files_src = glob('../src_snapshot/*.h5')
from msmbuilder.utils import map_drawn_samples
frames_abl = map_drawn_samples(abl_samples,traj_files_abl)
frames_src = map_drawn_samples(abl_samples,traj_files_abl)

for i in range(len(frames_abl)):
    msm_weight = msm_abl.populations_[i]
    frames_abl[i].save_pdb('abl_cluster_{0}_msm_weight_{1:.3f}.pdb'.format(i,msm_weight))

for i in range(len(frames_src)):
    msm_weight = msm_src.populations_[i]
    frames_abl[i].save_pdb('src_cluster_{0}_msm_weight_{1:.3f}.pdb'.format(i,msm_weight))
