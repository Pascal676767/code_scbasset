#!/bin/bash

# Job name:
#SBATCH --job-name=c5_s100_1000_enhancer

# Project:
#SBATCH --account=ec31

# Wall time limit:
#SBATCH --time=50:00:00
#SBATCH --ntasks=6
#SBATCH --cpus-per-task=1
#SBATCH --mem=250G

# Set the output of the slurmfile
#SBATCH --output=/fp/projects01/ec31/mathelier-group/pascalo/scBasset_internship/results_scbasset/log/output-%j.txt

## Set up job environment:
set -o errexit  # Exit the script on any error
set -o nounset  # Treat any unset variables as an error

module --quiet purge  # Reset the modules to the system default
module load Python/3.11.3-GCCcore-12.3.0
module list

##load the virtual environment
source /fp/projects01/ec31/mathelier-group/pascalo/scBasset_internship/Basset_tool/bin/activate

##go to the directory
cd /fp/projects01/ec31/mathelier-group/pascalo/scBasset_internship/ex_shap



## Do some work:
python ../code_scbasset/source/DeepSHAP/shap_multiprocessing.py data/enhancer/test_seqs.h5 data/enhancer/best_model.h5 5 1000 batch_enhancer_shuffle100_v2
