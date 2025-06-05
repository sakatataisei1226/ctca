#!/bin/bash
#PJM -L  "node=3"
#PJM --mpi proc=128
#PJM -L  "rscgrp=small"
#PJM -L  "elapse=10:00:00"
#PJM -x PJM_LLIO_GFSCACHE=/vol0004
#PJM -g [GROUP NAME]

. /vol0004/apps/oss/spack/share/spack/setup-env.sh
spack load /yhazdvl
spack load /upvlzyl

#Known issue: Path of dynamic link libraries of the operating system
export LD_LIBRARY_PATH=/lib64:$LD_LIBRARY_PATH

export EMSES_DEBUG=no

date

mpiexec ./mpiemses3D plasma.inp

date

# Postprocessing(visualization code, etc.)
