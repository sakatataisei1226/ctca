#!/bin/bash
#SBATCH -p gr20001a
#SBATCH -t 10:00
#SBATCH --rsc p=13:t=1:c=1
#SBATCH -o %x.%j.out

module load intel/2022.3

unset I_MPI_PMI_LIBRARY

mpiexec.hydra -n 4 ./sample3_requester : -n 1 ./sample3_coupler : -n 8 ./sample3_worker
date
