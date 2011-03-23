#!/bin/sh

# empty the folder
#make clean
#make clrout
#echo "Make clean"

PACKAGE_NAME=pkg.tgz
# tar all files necessary for compiling and running the program into a single file
tar cpzf $PACKAGE_NAME --exclude=result --exclude=trace --exclude=pkg.tgz --exclude=tags --exclude='*.o' --exclude=out --exclude='*.out' --exclude='*.swp' --exclude='.git' --exclude=submit.sh ./

# transfer the archive
Node1_dir=nical516:~/qsub
Node2_dir=chenwx-desktop:~/qsub
Node3_dir=~/qsub
scp $PACKAGE_NAME $Node1_dir
scp $PACKAGE_NAME $Node2_dir
scp $PACKAGE_NAME $Node3_dir
echo "Complete file transfer"

tar xpfz $PACKAGE_NAME -C ./

# submit jobs
qsub ccvil.pbs
