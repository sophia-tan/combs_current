#!/bin/bash

#export TMPDIR=/scratch
#export MYTMP=`mktemp -d`
#cd $MYTMP

#$ -S /bin/bash
#$ -o /netapp/home/skt/logs/integrin/
#$ -e /netapp/home/skt/logs/integrin/
#$ -pe smp 9
#$ -R y 
#$ -j y
#$ -l mem_free=18G
#$ -l arch=linux-x64
#$ -l h_rt=180:00:00
##$ -l netapp=2G

# Anything under here can be a bash script

date
hostname

cd ./Integrin_part_g/
echo "Working directory is $PWD"
cat submit.sh
python2 run.py
qstat -j $JOB_ID                                  # This is useful for debugging and usage purposes,
                                                  # e.g. "did my job exceed its memory request?"