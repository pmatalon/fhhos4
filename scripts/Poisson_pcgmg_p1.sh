#!/bin/bash

d=2
p=1

while getopts "d:p:" opt
do
   case "$opt" in
      d ) d="$OPTARG" ;;
      p ) p="$OPTARG" ;;
   esac
done

if [ $d = 2 ]
then
	listN=(32 64 128 256 512 1024)
elif [ $d = 3 ]
then
	if [ $p -le 2 ]
	then
		listN=(16 32 64 128)
	else
		listN=(8 16 32 64)
	fi
fi

for N in "${listN[@]}"
do
	echo "Poisson${d}D_pcgmg_n${N}_p${p}..."
	./bin/dghho -d $d -a s -s pcgmg -n $N -p $p > log/Poisson${d}D_pcgmg_n${N}_p${p}.txt
done
echo Finished.