#!/bin/bash

#SBATCH --account=paul-lab
#SBATCH --partition=neuro-largemem
#SBATCH --job-name=Criticality_Test
#SBATCH --output=crit_test_output.txt
#SBATCH --qos=medium
#SBATCH --ntasks=1
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=hgerm@brandeis.edu

#Load modules needed for the job
module load share_modules/MATLAB/R2018b

#Run your code
matlab --nodesktop -nosplash -nojvm -singleCompThread < criticality_tests_hpcc.m