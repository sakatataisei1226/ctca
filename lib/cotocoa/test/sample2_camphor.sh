#!/bin/bash
#SBATCH -p gr20001a
#SBATCH -t 10:00
#SBATCH --rsc p=29:t=1:c=1
#SBATCH -o %x.%j.out

module load intel/2022.3

unset I_MPI_PMI_LIBRARY

mpiexec.hydra -n 4 ./sample2_requester : -n 1 ./sample2_coupler : -n 8 ./sample2_worker1 : -n 8 ./sample2_worker2 : -n 8 ./sample2_worker3
date
