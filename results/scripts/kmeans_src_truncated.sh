#!/bin/bash
#
# walltime : maximum wall clock time (hh:mm:ss)
#PBS -l walltime=24:00:00
#
# join stdout and stderr
#PBS -j oe
#
# spool output immediately
#PBS -k oe
#
# specify queue
#PBS -q batch
#
# nodes: number of nodes
#   ppn: number of processes per node
#PBS -l nodes=1:ppn=8
#PBS -l mem=32G
# export all my environment variables to the job
#PBS -V
#
#
#
# filename for standard output (default = <job_name>.o<job_id>)
# at end of job, it is in directory from which qsub was executed
# remove extra ## from the line below if you want to name your own file
##PBS -o myoutput

# Change to working directory used for job submission
cd $PBS_O_WORKDIR

PATH=$PATH:'/cbio/jclab/home/fassj/anaconda2/bin'

source activate '/cbio/jclab/home/fassj/anaconda2'

# Launch my program.
python kmeans_truncated.py '../results/src_11401_tica.npz' 'src_11401_2000_t' 2000
