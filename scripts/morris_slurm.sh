#!/bin/bash
#SBATCH --job-name="sensitivity"
#SBATCH --partition=platinum
#SBATCH --time=200:00:00
#SBATCH --ntasks=50
#SBATCH --output=%j-%x-stdout.txt
#SBATCH --error=%j-%x-stderr.txt
#SBATCH --account=csd786
#SBATCH --qos=hcp-csd765
#SBATCH --nodes=1
#SBATCH --mem=200G

module load shared cpu/0.17.3
module load matlab/2022b-6eej3gx
cd /tscc/nfs/home/eafrancis
matlab -nodisplay -r "addpath(genpath('gitrepos/MATLAB')); SkelMuscleCa_Morris; quit"