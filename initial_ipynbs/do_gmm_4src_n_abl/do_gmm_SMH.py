"""
Build a GMM-MSM on tICA coordinates and plot the first two TICS with labels.
"""
import matplotlib
matplotlib.use('Agg')

from msmbuilder import example_datasets, cluster, msm, featurizer, lumping, utils, dataset, decomposition
from sklearn.pipeline import make_pipeline
import numpy as np
import matplotlib.pyplot as plt

tica_lagtime = 1600

dih = dataset.NumpyDirDataset("./dihedrals/")
X = dataset.dataset("./tica%d.h5" % tica_lagtime)
Xf = np.concatenate(X)

tica_model = utils.load("./tica%d.pkl" % tica_lagtime)
dih_model = utils.load("./dihedrals/model.pkl")

n_first = 2
n_components = 4

slicer = featurizer.FirstSlicer(n_first)
clusterer = cluster.GMM(n_components=n_components)
msm_model = msm.MarkovStateModel()

pipeline = make_pipeline(slicer, clusterer, msm_model)
s = pipeline.fit_transform(X)

p0 = make_pipeline(dih_model, tica_model, slicer)

trajectories = dataset.MDTrajDataset("./trajectories/*.h5")
selected_pairs_by_state = msm_model.draw_samples(s, 5)
samples = utils.map_drawn_samples(selected_pairs_by_state, trajectories)

for k, t in enumerate(samples):
    t.save("./pdbs/state%d.pdb" % k)


# Find the traj, frame indices for cluster center exemplars.
scores = map(lambda x: clusterer.score_samples(x[:, 0:n_first])[1], X)
max_scores = np.array(map(lambda x: x.max(0), scores))
traj_indices = max_scores.argmax(0)
frame_indices = [scores[trj_ind][:, i].argmax() for i, trj_ind in enumerate(traj_indices)]

selected_pairs_by_state = np.array(zip(traj_indices, frame_indices))[:, None]
samples = utils.map_drawn_samples(selected_pairs_by_state, trajectories)

for k, t in enumerate(samples):
    t.save("./pdbs/center%d.pdb" % k)

Y = p0.transform(samples)

plt.hexbin(Xf[:, 0], Xf[:, 1], bins='log')
plt.plot(clusterer.means_[:, 0], clusterer.means_[:, 1], 'k+', markersize=12, markeredgewidth=3)
map(lambda k: plt.annotate(k, xy=clusterer.means_[k, 0:2], fontsize=24), np.arange(n_components))
map(lambda y: plt.plot(y[:, 0], y[:, 1], 'x', markersize=8, markeredgewidth=2), Y)

plt.savefig('fig.png')
