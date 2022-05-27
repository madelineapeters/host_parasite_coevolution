#!/bin/bash
i=$1
g++ -std=gnu++11 continuousTimeMut_bash.cpp -o continuousTimeMut_bash

for r in {1..1300}
do
    ./continuousTimeMut_bash $i $r

if [ $r -eq 1 ]
then
	echo "r = 1"
	head -n 1 ./Data_${i}/data_${i}_${r}.csv > ./Data_${i}/data_final_${i}.csv
	head -n 1 ./Data_${i}/dataN_${i}_${r}.csv > ./Data_${i}/dataN_final_${i}.csv
fi

     tail -n 1 ./Data_${i}/data_${i}_${r}.csv >> ./Data_${i}/data_final_${i}.csv
     tail -n 1 ./Data_${i}/dataN_${i}_${r}.csv >> ./Data_${i}/dataN_final_${i}.csv
done

#head -n 1001 ./Data_${i}/data_final_${i}.csv > ./Data_${i}/data_final_${i}.csv
#head -n 1001 ./Data_${i}/dataN_final_${i}.csv > ./Data_${i}/dataN_final_${i}.csv

