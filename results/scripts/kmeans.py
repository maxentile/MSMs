# discretize and analyze!
import numpy as np
from msmbuilder import cluster
import cPickle
import pyemma

def discretize(path,name,n_clusters=500):
   ''' path is a string to an .npz file'''

   Y = [y[1] for y in np.load(path).items()]

   # run kmeeans clustering
   kmeans = pyemma.coordinates.cluster_mini_batch_kmeans(Y,k=n_clusters,max_iter=1000)
   dtrajs = [np.array(dtraj)[:,0] for dtraj in kmeans.get_output()]

#   kmeds = cluster.MiniBatchKMedoids(n_clusters=n_clusters,batch_size=n_clusters)
#   kmeds.fit([y[::downsampling_ratio] for y in Y])
#   dtrajs = kmeds.transform(Y)

   np.savez_compressed('{0}_dtrajs.npz'.format(name),*dtrajs)
   f = open('{0}_kmeds.pickle'.format(name),'w')
   cPickle.dump(kmeds,f)
   f.close()

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

   dtrajs = discretize(path=path,name=name,n_clusters=n_clusters)
