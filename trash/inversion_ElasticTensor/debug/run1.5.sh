#!/bin/bash
# Lines starting with #PBS are treated by bash as comments, but interpreted by qsub
# as arguments.  For more details about usage of these arguments see "man qsub"

#
# Name the job.

#PBS -N inv_vpvsCM1.4

#
# Set a walltime for the job. The time format is HH:MM:SS

# Run for 30 seconds:
#PBS -l walltime=1:00:00

#
#### running directory 
#PBS -d  /projects/jixi7887/work/US/inv_ET_BS

# Select one node, and only 1 processors per node
#you can see what properties are available
# with the "pbsnodes -a" command which will list all nodes and their properties

#PBS -l nodes=1:ppn=12

# # Join the Output and Errors in one file. you can also set the # path for this output.
# (see "man qsub" for details.)

#PBS -j oe

#
# Source the Dotkit init so you can pull in extra software via the "reuse" command:

. /curc/tools/utils/dkinit

#
# cd to the jobs working directory, which you can set above with a #PBS # directive,
# (see "man qsub" for details)

cd $PBS_O_WORKDIR

#
# Execute the program.
test_updateKeachjump_changeEtaSpace1.1_parall_Apr25_vpvsgdb2
#csh do1.5.csh temp_point4.txt 100 v14.2 
#csh do_ET_BS_CstMat_get_prior.csh temp_point1.txt 20 v14.2
#csh do_ET_BS_CstMat_mulP_Feb15Maincode_newV2L_changeInputc3.csh temp_point1.txt 50
# csh do_ET_BS_CstMat_get_prior.csh temp_point_t4.txt 1 v13.5
# This script needs to be submitted via qsub to run on the cluster
