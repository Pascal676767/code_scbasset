#!/bin/bash

# Job name:
#SBATCH --job-name=Create_training_data

# Project:
#SBATCH --account=ec31

# Wall time limit:
#SBATCH --time=01:00:00

#SBATCH --mem=16G

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
cd /fp/projects01/ec31/mathelier-group/pascalo/scBasset_internship/code_scbasset/source/preprocessing




## Do some work:
python SCAFE_make_csv.py 
