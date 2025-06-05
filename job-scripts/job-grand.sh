#!/bin/sh
#------ pjsub option --------#
#PJM -L rscgrp=ea
#PJM -L node=1
#PJM --mpi proc=32
#PJM -L elapse=1:00:00
#PJM -g [GROUP NAME]
#PJM -j
#------- Program execution -------#

date

export EMSES_DEBUG=no
mpiexec.hydra -n 32 mpiemses3D ./plasma.inp

date

# Postprocessing(visualization code, etc.)
