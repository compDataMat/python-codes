#!/bin/usr/env bash

s=("87690" "90000" "105000" "120000" "135000" "150000" "165000" "180000" "195000" "210000" "225000" "240000" "255000" "270000" "285000" "300000") 

for n in ${s[@]}; 
do
    echo Running input file : $n
    filenamedump="dump.stress.$n"
    filenamecube="$n.cube"
    python generate-stress-grids.py $filenamedump $filenamecube &
done
