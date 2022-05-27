#!/bin/bash
i=$1
g++ -std=gnu++11 continuousTimeMut_bash.cpp -o continuousTimeMut_bash

for r in {1..1000}
do
    ./continuousTimeMut_bash $i $r
done
