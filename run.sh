#!/bin/bash
for ((i = 1;i<=20;i=$i+1)) 
do
	sed "s/^functionToRun [0-9]*$/functionToRun $i/" configure.ini > temp
	cat temp > configure.ini
	echo $i
	cat configure.ini | grep '^functionToRun [0-9]*$'
	./submit.sh
	echo sleep
	sleep 15
done

rm temp

echo 'All job have been submitted successfully'
