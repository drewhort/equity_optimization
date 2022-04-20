#!/bin/bash

#SBATCH --job-name=no_le

#SBATCH --output=no_le.out.%j

#SBATCH --error=no_le.err.%j

#SBATCH -N 1

#SBATCH -p math-alderaan

singularity exec /storage/singularity/pyscipopt.sif /home/hortondr/equitable_facility_location/kpgrocery.sh
