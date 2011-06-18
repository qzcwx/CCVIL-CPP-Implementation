#!/bin/bash
make

for ((i = 1;i<=20;i=$i+1)) 
do
	for ((j = 4; j<=9; j=$j+1))
	do
		#	sed "s/^functionToRun [0-9]*$/functionToRun $i/" configure.ini > temp
		#	cat temp > configure.ini
		echo funcID $i, learnStrat $j
		#	cat configure.ini | grep '^functionToRun [0-9]*$'
		# ./submit.sh $i
		qsub ccvil.pbs -v jobID=$i,learnStrat=$j
		sleep 2
	done
done

echo 'All job have been submitted successfully'
