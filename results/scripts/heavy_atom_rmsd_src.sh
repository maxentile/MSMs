#!/bin/bash
# walltime : maximum wall clock time (hh:mm:ss)
#PBS -l walltime=72:00:00
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
#PBS -l nodes=1:ppn=16
#PBS -l mem=16G
#
# export all my environment variables to the job
#PBS -V
#
# job name (default = name of script file)
#PBS -N rmsd_abl
#
# mail settings (one or more characters)
# email is sent to local user, unless another email address is specified with PBS -M option 
# n: do not send mail
# a: send mail if job is aborted
# b: send mail when job begins execution
# e: send mail when job terminates
#PBS -m e
#
# filename for standard output (default = <job_name>.o<job_id>)
# at end of job, it is in directory from which qsub was executed
# remove extra ## from the line below if you want to name your own file
##PBS -o myoutput


cd $PBS_O_WORKDIR

PATH=$PATH:'/cbio/jclab/home/fassj/anaconda2/bin'

source activate '/cbio/jclab/home/fassj/anaconda2'

# Launch my program.
python ha_rmsd.py '../src_snapshot/*.h5' 'src_11401_10k_harmsd' 10000
