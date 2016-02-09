#!/bin/bash
#
# walltime : maximum wall clock time (hh:mm:ss)
#PBS -l walltime=0:30:00
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
#PBS -l nodes=1:ppn=1
#PBS -l mem=8G
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
python sanity_check.py

