#! /bin/bash

for r in {18..44}
do

arr_record1=( $(sed -n ${r}p inputMut_para.csv | cut -d ',' -f1) )
arr_record2=( $(sed -n ${r}p inputMut_para.csv | cut -d ',' -f2) )
arr_record3=( $(sed -n ${r}p inputMut_para.csv | cut -d ',' -f3) )
arr_record4=( $(sed -n ${r}p inputMut_para.csv | cut -d ',' -f4) )
arr_record5=( $(sed -n ${r}p inputMut_para.csv | cut -d ',' -f5) )
arr_record6=( $(sed -n ${r}p inputMut_para.csv | cut -d ',' -f6) )
arr_record7=( $(sed -n ${r}p inputMut_para.csv | cut -d ',' -f7) )
arr_record8=( $(sed -n ${r}p inputMut_para.csv | cut -d ',' -f8) )
arr_record9=( $(sed -n ${r}p inputMut_para.csv | cut -d ',' -f9) )
arr_record10=( $(sed -n ${r}p inputMut_para.csv | cut -d ',' -f10) )
arr_record11=( $(sed -n ${r}p inputMut_para.csv | cut -d ',' -f11) )
arr_record12=( $(sed -n ${r}p inputMut_para.csv | cut -d ',' -f12) )

bash inputMut_gen.sh ${arr_record1[${@}]} ${arr_record2[${@}]} ${arr_record3[${@}]} ${arr_record4[${@}]} ${arr_record5[${@}]} ${arr_record6[${@}]} ${arr_record7[${@}]} ${arr_record8[${@}]} ${arr_record9[${@}]} ${arr_record10[${@}]} ${arr_record11[${@}]} ${arr_record12[${@}]} ${r}
done
