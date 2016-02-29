# select clustering parameters
path='../abl_snapshot/*.h5'
n_clusters=10000
downsampling_ratio=5

# retrieve trajectory filenames
from glob import glob
filenames = glob(path)

# print mdtraj and msmbuilder versions
import mdtraj as md
print(md.version.full_version)
import msmbuilder.version
print(msmbuilder.version.full_version)

# select heavy atom indices
y = md.load_frame(filenames[0],index=0)
top = y.topology
heavy_atoms = top.select('mass>12')

# create MDTrajDataset
from msmbuilder import dataset
Y = dataset.MDTrajDataset(path=path,stride=downsampling_ratio,atom_indices=heavy_atoms)

# cluster
from msmbuilder import cluster
kmeds = cluster.MiniBatchKMedoids(n_clusters=n_clusters,batch_size=n_clusters,metric='rmsd')
dtrajs = kmeds.fit_transform(Y)

import numpy as np
