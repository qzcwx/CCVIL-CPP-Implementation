#!/bin/bash
make

# for ((i = 20;i>=1;i=$i-1)) 
for ((i=14;i<=18;i=$i+1)) 
do
	for ((j=0; j<=100; j=$j+10))
	do
		for ((k=2; k<=3; k=$k+1))
		do
			#	sed "s/^functionToRun [0-9]*$/functionToRun $i/" configure.ini > temp
			#	cat temp > configure.ini
			echo funcID $i,partialGroup $j,learnStrat $k
			#	cat configure.ini | grep '^functionToRun [0-9]*$'
			# ./submit.sh $i
			qsub ccvil.pbs -v jobID=$i,parGroup=$j,learnStrat=$k
			sleep 5
		done
	done
done

echo 'All job have been submitted successfully'
