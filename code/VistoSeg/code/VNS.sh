#!/bin/bash
#$ -cwd
#$ -l mem_free=100G,h_vmem=100G,h_stack=256M,h_fsize=100G
#$ -o /dcs04/lieber/marmaypag/lifespanDG_LIBD001/spatial_DG_lifespan/code/VistoSeg/code/logs/$TASK_ID_VNS.txt
#$ -e /dcs04/lieber/marmaypag/lifespanDG_LIBD001/spatial_DG_lifespan/code/VistoSeg/code/logs/$TASK_ID_VNS.txt
#$ -m e
#$ -M adramn83@gmail.com
#$ -t 1
#$ -tc 2

echo "**** Job starts ****"
date


echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"
echo "****"
echo "Sample id: V12N09-083_11_29_22_HRD_C1"
echo "****"

module load matlab/R2019a

toolbox='/dcs04/lieber/marmaypag/lifespanDG_LIBD001/spatial_DG_lifespan/code/VistoSeg/code'
fname='/dcs04/lieber/marmaypag/lifespanDG_LIBD001/spatial_DG_lifespan/processed-data/Images/VistoSeg/Capture_areas/V12N09-083_11_29_22_HRD_C1.tif'


matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), VNS('$fname',5)"

echo "**** Job ends ****"
date



