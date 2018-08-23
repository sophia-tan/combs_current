#!/bin/bash

#export TMPDIR=/scratch
#export MYTMP=`mktemp -d`
#cd $MYTMP

#$ -S /bin/bash
#$ -o /netapp/home/skt/logs/integrin/
#$ -e /netapp/home/skt/logs/integrin/
#$ -cwd
#$ -pe smp 10
#$ -j y
#$ -l mem_free=6G
#$ -l arch=linux-x64
##$ -l netapp=2G

# Anything under here can be a bash script

date
hostname

python2 run.py
qstat -j $JOB_ID                                  # This is useful for debugging and usage purposes,
                                                  # e.g. "did my job exceed its memory request?"
