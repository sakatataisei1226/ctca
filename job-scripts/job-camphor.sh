#!/bin/bash
#SBATCH -p [QUEUE NAME]
#SBATCH --rsc p=32:t=1:c=1
#SBATCH -t 1:00:00
#SBATCH -o stdout.%J.log
#SBATCH -e stderr.%J.log

set -x

module load intel/2022.3 intelmpi/2022.3
module load hdf5/1.12.2_intel-2022.3-impi
module load fftw/3.3.10_intel-2022.3-impi

export EMSES_DEBUG=no

date

srun ./mpiemses3D plasma.inp

date

# Postprocessing(visualization code, etc.)
