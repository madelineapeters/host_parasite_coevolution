#!/bin/bash

x=$1

g++ -std=gnu++11 continuousTimeIsland_bash.cpp -o continuousTimeIsland_bash

for r in {1..46}
do
    ./continuousTimeIsland_bash $x $r
done
