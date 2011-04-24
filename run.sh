#!/bin/bash
make

for ((i = 20;i>=1;i=$i-1)) 
#for ((i=14;i<=18;i=$i+1)) 
do
			#	sed "s/^functionToRun [0-9]*$/functionToRun $i/" configure.ini > temp
			#	cat temp > configure.ini
			echo funcID $i
			#	cat configure.ini | grep '^functionToRun [0-9]*$'
			# ./submit.sh $i
			qsub ccvil.pbs -v jobID=$i
		#	sleep 1
done

echo 'All job have been submitted successfully'
