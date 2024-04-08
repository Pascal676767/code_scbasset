#!/bin/bash

# Job name:
#SBATCH --job-name=Train_HP_notfilt

# Project:
#SBATCH --account=ec31

# Wall time limit:
#SBATCH --time=09:00:00

#SBATCH --mem=8G
#SBATCH --partition=accel --gres=gpu:1

# Set the output of the slurmfile
#SBATCH --output=/fp/projects01/ec31/mathelier-group/pascalo/scBasset_internship/results_scbasset/log/output-%j.txt

## Set up job environment:
set -o errexit  # Exit the script on any error
set -o nounset  # Treat any unset variables as an error

module --quiet purge  # Reset the modules to the system default
module load Python/3.11.3-GCCcore-12.3.0
## module load tensorflow-probability/0.19.0-foss-2022a-CUDA-11.7.0
## module load h5py/3.6.0-foss-2021b
module list

##load the virtual environment
source /fp/projects01/ec31/mathelier-group/pascalo/scBasset_internship/Basset_tool/bin/activate

##go to the directory
cd /fp/projects01/ec31/mathelier-group/pascalo/scBasset_internship/results_scbasset/enhancer_final/processed/



## Do some work:

python /fp/projects01/ec31/mathelier-group/pascalo/scBasset_internship/code_scbasset/source/scbasset_train.py --input_folder filtered/  --epochs  1000 --out_path ../../../results_scbasset/enhancer_final/model/l2_regularization/double/filtered
