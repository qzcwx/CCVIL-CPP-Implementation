#!/bin/bash
## This is a script written for submitting a job to cluster

# empty the folder
#make clean
#make clrout
#echo "Make clean"

# tar all files necessary for compiling and running the program into a single file
#PACKAGE_NAME=pkg.tgz
tar cpzf pkg.tgz --exclude=result --exclude=trace --exclude=pkg.tgz --exclude=tags --exclude='*.o' --exclude=out --exclude='*.out' --exclude='*.swp' --exclude='.git' --exclude='Run-CCVIL.*' --exclude=submit.sh --exclude=run.sh ./
#tar cpzf $PACKAGE_NAME --exclude=result --exclude=trace --exclude=pkg.tgz --exclude=tags --exclude='*.o' --exclude=out --exclude='*.out' --exclude='*.swp' --exclude='.git' --exclude=submit.sh ./

# transfer the archive
Node1_dir=nical516:~/qsub
Node2_dir=chenwx-desktop:~/qsub
#Node2_dir=~/qsub
scp pkg.tgz $Node1_dir
scp pkg.tgz $Node2_dir
#scp $PACKAGE_NAME $Node3_dir

echo "Complete file transfer"

echo Job ID $1

# submit jobs
qsub ccvil.pbs -v job_id=$1
