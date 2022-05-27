#!/bin/bash
i=$1
g++ -std=gnu++11 continuousTimeIsland_bash.cpp -o continuousTimeIsland_bash

for r in {1..1000}
do
    ./continuousTimeIsland_bash $i $r

if [ $r -eq 1 ]
then
	echo "r = 1"
	head -n 1 ./Data_Island_${i}/data_island_${i}_${r}.csv > ./Data_Island_${i}/data_final_${i}.csv
	head -n 1 ./Data_Island_${i}/dataN_island_${i}_${r}.csv > ./Data_Island_${i}/dataN_final_${i}.csv
fi

     tail -n 1 ./Data_Island_${i}/data_island_${i}_${r}.csv >> ./Data_Island_${i}/data_final_${i}.csv
     tail -n 1 ./Data_Island_${i}/dataN_island_${i}_${r}.csv >> ./Data_Island_${i}/dataN_final_${i}.csv
done
