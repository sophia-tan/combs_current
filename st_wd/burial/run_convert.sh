#!/bin/bash

cat ifgs.txt | while read a
do
python gen_scores.py $a > "${a}log" & 
done
