#!/bin/bash
i=$1
g++ -std=gnu++11 continuousTimeMut_bash.cpp -o continuousTimeMut_bash

for r in {1..1300}
do

 awk -F '","'  'BEGIN {OFS=","} { if (toupper($1) <= 150 )  print }' ./Data_${i}/data_${i}_${r}.csv > ./Data_${i}/data_${i}_${r}.csv
awk -F '","'  'BEGIN {OFS=","} { if (toupper($1) <= 150 )  print }' ./Data_${i}/dataN_${i}_${r}.csv > ./Data_${i}/dataN_${i}_${r}.csv

if [ $r -eq 1 ]
then
	echo "r = 1"
	head -n 1 ./Data_${i}/data_${i}_${r}.csv > ./Data_${i}/data_final_${i}.csv
	head -n 1 ./Data_${i}/dataN_${i}_${r}.csv > ./Data_${i}/dataN_final_${i}.csv
fi

     tail -n 1 ./Data_${i}/data_${i}_${r}.csv >> ./Data_${i}/data_final_${i}.csv
     tail -n 1 ./Data_${i}/dataN_${i}_${r}.csv >> ./Data_${i}/dataN_final_${i}.csv
done

head -n 1001 ./Data_${i}/data_final_${i}.csv > ./Data_${i}/data_final_${i}.csv
head -n 1001 ./Data_${i}/dataN_final_${i}.csv > ./Data_${i}/dataN_final_${i}.csv

