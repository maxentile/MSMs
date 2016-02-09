import mdtraj as md
import pyemma
import msmbuilder
import numpy as np
from glob import glob
import cPickle

#import matplotlib.pyplot as plt
#plt.rc('font', family='serif')

def featurize(path_to_files,model_name):

    files = glob(path_to_files)
    print("Number of files: {0}".format(len(files)))

    # timestep between frames: 250 picoseconds

    #### A. FEATURE EXTRACTION ####

    ### Step 1: identifying interresidue contacts that change

    # compute full contact maps for a strided subset of the simulation frames

    strided_distances=[]
    stride=100 # stride within trajectory
    traj_thin=5 # only look at 1 in traj_thin trajectories
    scheme = 'ca'
    threshold = 0.8 # contact threshold in angstroms

    for f in files[::traj_thin]:
        traj = md.load(f,stride=stride)
        distances,residue_pairs = md.compute_contacts(traj,scheme=scheme)
        strided_distances.append(distances)

    strided_distances = np.vstack(strided_distances)

    # identify contacts that change by counting how many times the distances were
    # greater than and less than the threshold
    num_times_greater_than = (strided_distances>threshold).sum(0)
    num_times_less_than = (strided_distances<threshold).sum(0)
    changed = (num_times_greater_than > 0) * (num_times_less_than > 0)
    print("Number of contacts that changed: {0}".format(changed.sum()))
    print("Total number of possible contacts: {0}".format(len(residue_pairs)))

    # now turn this bitmask into a list of relevant residue pairs
    respairs_that_changed = residue_pairs[changed]

    # save this list!
    np.save('{0}_respairs_that_changed.npy'.format(model_name),respairs_that_changed)

    ### Step 2: extract these selected features from the full dataset

    X = []

    traj_thin=1 # only look at 1 in traj_thin trajectories
    files_of_interest = files[::traj_thin]

    for i,f in enumerate(files_of_interest):
        print('{0}/{1}'.format(i,len(files_of_interest)))
        traj = md.load(f)
        distances,_ = md.compute_contacts(traj,contacts=respairs_that_changed,scheme=scheme)

        X.append(distances)

    print("Initial dimensionality: {0}".format(X[0].shape[1]))
    print("# frames: {0}".format(np.vstack(X).shape[0]))

    ##### B. KINETIC DISTANCE LEARNING #####
    tica = pyemma.coordinates.tica(X)
    Y = tica.get_output()

    # save tica model and output
    np.savez_compressed('{0}_tica.npz'.format(model_name),*Y)
    print("Dimensionality after tICA, retaining enough eigenvectors to explain 0.95 of kinetic variation: {0}".format(np.vstack(Y).shape[1]))


if __name__=='__main__':
    import sys
    path = sys.argv[1]
    model_name=sys.argv[2]

    featurize(path,model_name)
