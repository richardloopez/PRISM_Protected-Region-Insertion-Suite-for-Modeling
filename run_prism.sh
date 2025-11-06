#!/bin/bash
#SBATCH --job-name=PRISM_MODELLER_JOB
#SBATCH --output=slurm_output_%j.out
#SBATCH --error=slurm_error_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16

# Load required modules
module load python/3.8.12

# This variable will be read by config.py
export SLURM_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK

# Do NOT cd into the 'prism/' package directory
# Run the module from the project root
echo "Starting PRISM Pipeline from $PWD..."
srun python3 -m prism.controller > prism_pipeline.log
echo "PRISM Pipeline finished. Check prism_pipeline.log for details."
