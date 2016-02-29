# discretize and analyze!
import numpy as np
from msmbuilder import cluster
import cPickle
import mdtraj as md

def discretize(path,name,n_clusters=500,downsampling_ratio=5):
   from time import time
   t = time()
   from glob import glob
   filenames = glob(path)

   # load heavy-atom coordinates from all trajectories
   y = md.load_frame(filenames[0],index=0)
   top = y.topology
   heavy_atoms = top.select('mass>=12')
#   Y = [y.atom_slice(heavy_atoms) for y in md.load(filenames,stride=downsampling_ratio)]
   from msmbuilder import dataset
   Y = dataset.MDTrajDataset(path=path,stride=downsampling_ratio,atom_indices=heavy_atoms)

   print('Loaded trajectories!')
   print(time() - t)

   # run kmedoids clustering on downsampled data
   kmeds = cluster.MiniBatchKMedoids(n_clusters=n_clusters,max_iter=1000,batch_size=n_clusters,metric='rmsd')
   dtrajs = kmeds.fit_transform(Y)
   print('Finished mini-batch k-medoids!')
   print(time() - t)

#   # clear Y from memory
#   del(Y)
#   print('Cleared trajectories!')

#   # assign all frames to nearest cluster generators
#   dtrajs = []
#   for filename in filenames:
#      y = md.load(filename)
#      dtrajs.append(kmeds.transform([y.atom_slice(heavy_atoms)])[0])
#   del(y)
#   print('Discretized trajectories!')
#   print(time() - t)

   # save clustering results
   np.savez_compressed('{0}_dtrajs.npz'.format(name),*dtrajs)
   print('Saved discrete trajectories!')
   print(time() - t)

   # save clustering model
   f = open('{0}_kmeds.pickle'.format(name),'w')
   cPickle.dump(kmeds,f)
   f.close()
   print('Pickled k-medoids model!')
   print(time() - t)

   return dtrajs

if __name__=='__main__':
   import sys
   path=sys.argv[1]
   name=sys.argv[2]
   try:
      n_clusters=int(sys.argv[3])
   except:
      n_clusters = 500

   print(path)
   print(name)
   dtrajs = discretize(path,name=name,n_clusters=n_clusters,downsampling_ratio=1)
