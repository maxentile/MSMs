There are a lot of random scripts here, corresponding to stages of MSM building I ran separately and inspected before proceeding.

The shell scripts are all just one-liners*, calling a single python script with a set of arguments, plus requesting the appropriate cluster resources.

The python scripts include several attempts at discretization, and some assorted analysis stuff.
- featurize.py -- Converted each frame into a feature vector of inter-residue distances (only looking at the subset of these distances that crossed some threshold over the course of the simulation), then projected onto the leading tICs
- discretize.py -- First pass at discretizing the resulting tICA projection -- used minibatch kmedoids, but had an insufficient maximum number of iterations
- kmeans.py -- Second pass at discretizing the tICA projection above -- used minibatch kmeans, but with a reasonable iteration budget
- kmeans_truncated.py -- Third pass at discretizing the tICA projection above -- to reduce compuational cost / avoid timing out, use only the top-50 tICs.
- ha_rmsd.py -- Heavy-atom RMSD Minibatch K-Medoids -- ignored the tICA projection, and performed clustering using heavy-atom RMSD instead
- analyze_dtrajs.py and analyze_dtrajs_no_plotting.py -- Various analyses, including identifying the number of metastable states, plotting implied timescales, and a couple sanity checks
