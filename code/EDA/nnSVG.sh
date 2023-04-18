#!/bin/bash
#$ -cwd
#$ -l mem_free=15G,h_vmem=15G,h_fsize=100G
#$ -pe local 16
#$ -N nnSVG
#$ -o logs/nnSVG.txt
#$ -e logs/nnSVG.txt
#$ -m e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"

## Load the R module
module load conda_R/4.2

## List current modules for reproducibility
module list

## Edit with your job command
Rscript nnSVG.R

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/