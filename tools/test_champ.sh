#!/bin/bash
#SBATCH -p short                   # partition (queue)
##SBATCH -p normal                   # partition (queue)
##SBATCH -N 2                      # number of nodes
#SBATCH --exclude=tcn1,tcn2,tcn5
#SBATCH -n 5                     # number of cores
##SBATCH --mem 20000                 # memory pool for all cores
#SBATCH -t 00:05:00                 # time (D-HH:MM)
#SBATCH -o slurm.%N.%j.out        # STDOUT
#SBATCH -e slurm.%N.%j.err        # STDERR
##SBATCH --mail-user=p.l.t.lopeztarifa@vu.nl # send-to address
##SBATCH --dependency=afterok:57135
#
module load pre2019
module load intel/2018b mpi/impi/18.0.4
SRC=/home/plopez/Programs/champ/bin
JOBNAME=$(echo $1 | sed s/\.inp//g)
echo "JOBNAME is: $JOBNAME" 
######## LOAD MODULES  #########
date
 ls -l
 beginning=`date +%s`
 srun  $SRC/vmc.mov1  < $JOBNAME.inp > $JOBNAME.out 
 end=`date +%s`
 ls -l
date
#-------------------------------TIME----------------------------------------------#
  time=`expr $end "-" $beginning`
  hour=`echo "($time)/60/60"|bc`
  min=`echo "($time-$hour*60*60)/60"|bc`
  sec=`echo "$time-$min*60-$hour*60*60"|bc`
  echo "---------------------------------------" >> $JOBNAME.out
  echo "Job execution "human" time:" >> $JOBNAME.out
  echo $hour hours $min minutes $sec seconds >> $JOBNAME.out
  echo "---------------------------------------" >> $JOBNAME.out
#
