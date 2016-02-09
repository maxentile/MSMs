# reads in discrete trajectories, generates a report and writes out
# the following figures:
# - histogram of trajectory lengths
# - implied timescales

# to-do:
# - eigenvalue spectrum
# - cktest

import numpy as np
import numpy.random as npr
import pyemma.msm as msm
import pyemma.plots as mplt
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('Agg')
plt.rc('font',family='serif')

def load_dtrajs(filename):
   ''' Loads discrete trajectories from ``.npz`` files.

   Thinly wraps ``numpy.load``, discarding array names.

   Parameters
   ----------
   filename : string
     The file to read.

   Returns
   -------
   dtrajs : list of arrays
     Each array in the list is flat and contains integers.
   '''
   dtrajs = [y[1] for y in np.load(filename).items()]
   return dtrajs

def plot_traj_length_distribution(dtrajs,bins=20):
   ''' Create a histogram of trajectory lengths. '''
   plt.hist([len(x) for x in dtrajs],bins=20);
   plt.xlabel('Trajectory length (snapshots)')
   plt.ylabel('Number ')

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


def plot_implied_timescales(dtrajs,nits=15,model_name=''):
   ''' Compute and plot implied timescales.

   Parameters
   ----------
   dtrajs
   nits
   model_name

   Returns
   -------
   its : :class:`ImpliedTimescales <pyemma.msm.estimators.implied_timescales.ImpliedTimescales>` object
   '''

   its = msm.its(dtrajs, lags=lags, nits=nits,errors='bayes')
   mplt.plot_implied_timescales(its,dt=0.25,units='ns')
   plt.title(model_name)
   return its


def find_lag_time(its):
   ''' TO-DO: use implied timescales to select lag-time automatically'''
   pass

def plot_microstate_free_energies(msm,model_name,top_k=50):
   ''' Use stationary distribution of an MSM to compute the relative
   free energies of the top-k most stable microstates.

   Parameters
   ----------
   msm : pyemma.msm object
   model_name : string
   top_k : int
   '''
   f_i = -np.log(sorted(msm.stationary_distribution))[::-1];
   f_i -= f_i.min();

   plt.plot(f_i[:top_k], '.')
   plt.ylabel(r'Relative free energy ($k_B T \ln \pi_i$)')
   plt.xlabel('Microstate index ($i$)')
   plt.title(model_name)
   plt.savefig('{0}_microstate_free_energies.pdf'.format(model_name))
   plt.close()

def estimate_n_macrostates(msm,metastability_threshold):
   ''' Compute implied timescales at the given lag_time, return
   the number of processes with lower relative free energy
   than the metastability_threshold.

   This number of slow processes should roughly correspond to the number of
   metastable macrostates.

   Parameters
   ----------
   msm : pyemma.msm.MarkovStateModel object
     estimated
   metastability_threshold : float
     unitless (multiples of k_B T)

   Returns
   -------
   n_processes : int
   '''

   f_i = -np.log(sorted(msm.stationary_distribution))[::-1];
   f_i -= f_i.min();
   return np.sum(f_i<metastability_threshold)

def plot_spectrum(msm, model_name,dt,metastability_threshold=None):
   ''' create a simple plot of the eigenvalues of the msm transition matrix
   - represent each eigenvalue by a horizontal line
   - make vertical axis log-scale
   - label the vertical axis with both implied timescale (left) and relative free energy
   - (optionally draw the metastability threshold on this plot)
   '''
   pass


def find_most_metastable_clusters(M,top_k=10):
   '''
   Without performing any coarse-graining, identify the top-k most metastable
   microstates.


   Parameters
   ----------
   M : pyemma.msm object
      estimated MSM
   top_k : int
      number of microstates to return

   Returns
   -------
   clusters : list of tuples
      each element of clusters includes:
         cluster index (int)
         boltzmann weight (float)
         count fraction (float)
   '''
   flat_dtrajs = np.hstack(M.discrete_trajectories_full)
   cluster_inds = sorted(np.arange(len(M.stationary_distribution)),key=lambda i:-M.stationary_distribution[i])

   clusters = []
   for i in range(top_k):
      ind = cluster_inds[i]
      clusters.append((ind, # cluster index
           M.stationary_distribution[ind], # boltzmann weight
          1.0*sum(flat_dtrajs==ind)/len(flat_dtrajs))) # count fraction
   return clusters

def from_flattened_label_index_to_traj_frame_pair(label_index,trajs):
   ''' If you've taken a list of 1D arrays and concatentated them end-to-end
   to form a new array, and you know the index of an array element in that new
   flattened array, what is the corresponding (list_ind,element_ind) pair?
   '''
   len_trajs = [len(traj) for traj in trajs]
   cumsum = 0
   for i,traj in enumerate(trajs):
      if cumsum + len(traj) > label_index:
         return i,label_index-cumsum
      else:
         cumsum += len(traj)
   raise Exception('label_index is too large!')


def coarse_grain(dtrajs,n_metastable):
   ''' perform a coarse-graining using PCCA+'''

   from msmbuilder import lumping
   pcca = lumping.PCCAPlus(n_macrostates=n_metastable)
   pcca.fit(dtrajs)
   return pcca

def draw_exemplars(pcca_obj,n_exemplars_per_state=30):
   ''' Given a coarse-grained model, draw exemplar trajectory_id's and
   frame indices

   WARNING: this has been giving nonsensical output:
     - Some of the frame indices are larger than the length of their corresponding trajectories
     - Variation among exemplars within a single cluster is large relative to
      variation among exemplars from different clusters (e.g. containing both
      DFG-in and DFG-out conformations)
     - Best guess as to the cause of this behavior: the number / order of trajectory files
     in the munging directory are changing over time. If this script is called after
     featurizing/discretization, then the state of the directory could have changed
       - Possible workarounds:
         - Filter by edit date? (using the edit date of the )
         - When doing the feature-extraction step, record the filename in a dictionary
   '''


   pass

def analyze_dtrajs(filename,
               traj_files,
               model_name,
               lags,
               lag=100,
               dt=0.25,
               metastability_threshold=6):
   ''' Applies full analysis / coarse-graining pipeline. '''

   ### 1. SETUP / LOADING
   # initialize report file
   report_filename = '{0}_report.txt'.format(model_name)
   report = open(report_filename,'w')

   # load dtrajs from file
   raw_dtrajs = load_dtrajs(filename)

   # examine length distribution of these trajectories
   plot_traj_length_distribution(raw_dtrajs)

   # trim unusably short trajectories
   dtrajs = trim_dtrajs(raw_dtrajs)

   ### 2. ESTIMATE AND ANALYZE MSM
   # compute and plot implied timescales
   its = plot_implied_timescales(dtrajs,nits=30,model_name=model_name)

   # select a lag time (for now manual, soon automatic: see ``find_lag_time``
   #lag=find_lag_time(its)

   # estimate MSM
   M = msm.MaximumLikelihoodMSM(lag)
   M.fit(dtrajs)

   # estimate the number of metastable states
   n_metastable_states = estimate_n_macrostates(M,metastability_threshold)
   report.write('# metastable states at metastability threshold {0} kT : {1}\n'.format(
                metastability_threshold,n_metastable_states))

   ### 3. FIND AND SAVE EXEMPLAR STRUCTURES
   exemplars_per_cluster=30

   ### 3.a. Just use the most metastable microstates
   # 3.a.1. Identify metastable microstates and their relative populations
   top_k=20
   metastable_microstates = find_most_metastable_clusters(M,top_k=top_k)
   report.write('Top-{0} most metastable stable microstates: '.format(top_k))
   report.writelines(['  '+str(m) for m in metastable_microstates])

   # 3.a.2. Draw samples from each microstate and save them out as PDBs
   microstate_ids = [m[0] for m in metastable_microstates]
   dtrajs_flat = np.hstack(raw_dtrajs)
   for i in microstate_ids:
      flat_inds = np.where(dtrajs_flat==i)[0]
      npr.shuffle(flat_inds)
      select_inds = flat_inds[:exemplars_per_cluster]
      traj_frame_pairs = [from_flattened_label_index_to_traj_frame_pair(ind,dtrajs) for ind in select_inds]
      exemplars = [md.load_frame(traj_files[traj],frame) for (traj,frame) in traj_frame_pairs]
      exemplars_traj = exemplars[0].join(exemplars[1:])
      exemplars_traj.superpose(exemplars_traj)
      exemplars_traj.save_hdf5('{0}_cluster_{1}.h5'.format(model_name,i))

   #### 3.b. TO-DO: USE COARSE-GRAINED MODEL INSTEAD
   ## 3.b.1 coarse-grain to the estimated number of metastable states
   #pcca = coarse_grain(dtrajs,n_metastable=n_metastable_states)

   ## 3.b.2 draw exemplar structures from the coarse-grained model
   #draw_exemplars(pcca_obj)


   report.close()


if __name__=='__main__':
   # parse command line arguments:
   import sys
   from glob import glob

   # path to discrete trajectories
   try:
     filename = sys.argv[1]
   except:
     abl_data = '../results/src_11401_2000_dtrajs.npz'
     filename=abl_data

   # path to actual trajectories
   try:
     trajectories = glob(sys.argv[2])
   except:
     trajectories=glob('../src_snapshot/*.h5')

   # model name
   try:
     model_name = sys.argv[3]
   except:
     model_name='src'

   # lag-time for reversible estimation
   try:
      lag=int(sys.argv[4])
   except:
      lag=50

   # metastability_threshold (as a unitless relative free energy)
   try:
      metastability_threshold=float(sys.argv[5])
   except:
      metastability_threshold=6


   print('File name: ',filename)
   print('Model name: ', model_name)
   print('Lag: ',lag)
   print('Metastability threshold:', metastability_threshold)

   lags = [1,2,5,10,20,50,100,200,300,400,500]

   analyze_dtrajs(filename,trajectories,model_name,lags,lag=lag,dt=0.25,
                  metastability_threshold=metastability_threshold)
