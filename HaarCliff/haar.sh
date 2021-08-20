#!/bin/bash
#SBATCH --account=def-rgmelko
#SBATCH --mem=4G
#SBATCH --time=4:00:00
#SBATCH --mail-user=rokas@mit.edu
#SBATCH --mail-type=ALL

module add julia
julia ~/jstuff/haar_cliff_mpt/HaarCliff/haar1.jl
