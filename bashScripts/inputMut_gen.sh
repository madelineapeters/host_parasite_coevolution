#!/bin/bash
maxtime=$1
kappaH=$2
kappaP=$3
phi=$4
omega=$5
beta11=$6
beta12=$7
beta21=$7
beta22=$6
alpha=$8
delta=$9
gamma=${10}
mu=${11}
nu=${12}

i=${13}

echo "Input file for main.cpp" > inputMut_${i}.txt

VARtime="Maximum time (tmax): "
VARtime+="${maxtime}"

VARkappaH="Host population size (kappaH): "
VARkappaH+="${kappaH}"

VARkappaP="Parasite population size (kappaP): "
VARkappaP+="${kappaP}"

VARphi="Host migration rate (phi): "
VARphi+="${phi}"

VARomega="Parasite migration rate (omega): "
VARomega+="${omega}"

VARbeta="Transmission matrix (beta): "
VARbeta1="${beta11}"
VARbeta1+=" "
VARbeta1+="${beta12}"
VARbeta2="${beta21}"
VARbeta2+=" "
VARbeta2+="${beta22}"

VARalpha="Virulence (alpha): "
VARalpha+="${alpha}"

VARdelta="Host turnover rate (delta): "
VARdelta+="${delta}"

VARgamma="Parasite turnover rate (gamma): "
VARgamma+="${gamma}"

VARmu="Host mutation rate (mu): "
VARmu+="${mu}"

VARnu="Parasite mutation rate (nu): "
VARnu+="${nu}"

echo $VARtime >> inputMut_${i}.txt
echo $VARkappaH >> inputMut_${i}.txt
echo $VARkappaP >> inputMut_${i}.txt
echo $VARphi >> inputMut_${i}.txt
echo $VARomega >> inputMut_${i}.txt
echo $VARbeta >> inputMut_${i}.txt
echo $VARbeta1 >> inputMut_${i}.txt
echo $VARbeta2 >> inputMut_${i}.txt
echo $VARalpha >> inputMut_${i}.txt
echo $VARdelta >> inputMut_${i}.txt
echo $VARgamma >> inputMut_${i}.txt
echo $VARmu >> inputMut_${i}.txt
echo $VARnu >> inputMut_${i}.txt
