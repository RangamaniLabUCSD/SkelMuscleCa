#!/bin/bash
#SBATCH --job-name="mito_fit"
#SBATCH --output=%j-%x-stdout.txt
#SBATCH --error=%j-%x-stderr.txt
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=30
#SBATCH --mem=200G
#SBATCH --account=ddp494
#SBATCH --export=ALL
#SBATCH -t 48:00:00

module load cpu/0.17.3b
module load matlab/2022b/lefe4oq
cd /home/efrancis
matlab -nodisplay -r "addpath(genpath('gitrepos/MATLAB')); SkelMuscleCa_EstimateDriver; quit"