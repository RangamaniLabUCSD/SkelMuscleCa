#!/bin/bash
#SBATCH --job-name="mito_est"
#SBATCH --partition=condo
#SBATCH --time=200:00:00
#SBATCH --ntasks=30
#SBATCH --output=%j-%x-stdout.txt
#SBATCH --error=%j-%x-stderr.txt
#SBATCH --account=csd786
#SBATCH --qos=condo
#SBATCH --nodes=1
#SBATCH --mem=100G

module load shared cpu/0.17.3
module load matlab/2022b-6eej3gx
cd /tscc/nfs/home/eafrancis
matlab -nodisplay -r "addpath(genpath('gitrepos/MATLAB')); SkelMuscleCa_fitMitoParam; quit"