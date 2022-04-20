#!/bin/bash

#SBATCH --job-name=no_min

#SBATCH --output=no_min.out.%j

#SBATCH --error=no_min.err.%j

#SBATCH -N 1

#SBATCH -p math-alderaan

singularity exec /storage/singularity/pyscipopt.sif /home/hortondr/equitable_facility_location/kpgrocery.sh
