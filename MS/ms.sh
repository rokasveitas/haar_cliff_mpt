#!/bin/bash
#SBATCH --account=def-rgmelko
#SBATCH --mem=124G
#SBATCH --time=5-00:00:00
#SBATCH --mail-user=rokas@mit.edu
#SBATCH --mail-type=ALL
#SBATCH --array=1-384

module add julia
julia ~/jstuff/ms.jl $SLURM_ARRAY_TASK_ID
