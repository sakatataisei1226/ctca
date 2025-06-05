#!/bin/bash
#SBATCH -p gr20001a
#SBATCH -t 10:00
#SBATCH --rsc p=9:t=1:c=1
#SBATCH -o %x.%j.out

export PYTHONUNBUFFERED=1
set -x

pip install mpi4py

srun --multi-prog multi4.conf
date
