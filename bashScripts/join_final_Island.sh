#!/bin/bash
x=$1
for i in {45..53}
do

if [ $i -eq 18 ]
then
        echo "initializing final .csv"
	chmod 777 ./Data_Island_${i}/data_final_${i}.csv
	chmod 777 ./Data_Island_${i}/dataN_final_${i}.csv

        cat ./Data_Island_${i}/data_final_${i}.csv > ./data_final_${x}.csv
        cat ./Data_Island_${i}/dataN_final_${i}.csv > ./dataN_final_${x}.csv
	echo ${i}
	wc -l ./data_final_${x}.csv
else
	chmod 777 ./Data_Island_${i}/data_final_${i}.csv
        chmod 777 ./Data_Island_${i}/dataN_final_${i}.csv

        cat ./Data_Island_${i}/data_final_${i}.csv >> ./data_final_${x}.csv
        cat ./Data_Island_${i}/dataN_final_${i}.csv >> ./dataN_final_${x}.csv
	echo ${i}
        wc -l ./data_final_${x}.csv
fi
done
