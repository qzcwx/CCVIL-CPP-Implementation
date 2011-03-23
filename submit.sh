#!/bin/sh

# empty the folder
make clean
make clrout
echo "Make clean"

# transfer files
Node1_dir=nical516:~/qsub
Node2_dir=chenwx-desktop:~/qsub
Node3_dir=~/qsub
scp -r ../CCVIL $Node1_dir
scp -r ../CCVIL $Node2_dir
scp -r ../CCVIL $Node3_dir
echo "Complete file transfer"

# submit jobs
qsub ccvil.pbs
