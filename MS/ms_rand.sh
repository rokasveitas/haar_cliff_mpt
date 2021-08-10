#!/bin/bash
#SBATCH --account=def-rgmelko
#SBATCH --mem=120G
#SBATCH --time=1-00:00:00
#SBATCH --mail-user=rokas@mit.edu
#SBATCH --mail-type=ALL
#SBATCH --array=108-192

module add julia
julia ~/jstuff/haar_cliff_mpt/MS/ms_rand.jl $SLURM_ARRAY_TASK_ID
