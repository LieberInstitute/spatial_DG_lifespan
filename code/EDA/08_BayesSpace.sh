#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=200G,h_vmem=200G,h_fsize=200G
#$ -N bayesSpace_k_8
#$ -o logs/bayesSpace_k_8.txt
#$ -e logs/bayesSpace_k_8.txt
#$ -m e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"

## Load the R module (absent since the JHPCE upgrade to CentOS v7)
module load conda_R/4.2

## List current modules for reproducibility
module list

## Edit with your job command
Rscript 08_BayesSpace.R

echo "**** Job ends ****"
date