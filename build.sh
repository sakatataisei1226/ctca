#!/bin/bash

# Environment variable cotocoa_mode determines the mode of EMSES.
if [ -z "$cotocoa_mode" ]; then
	echo "Environment variable cotocoa_mode must be one of the following: {COTOCOA_DISABLED, COTOCOA_REQUESTER, COTOCOA_WORKER}. Exiting."
	exit 1
fi

hostname=$(uname -n)

if [[ $hostname == fn* ]]; then # on Fugaku
    hostname=Fugaku
	echo "The current environment has been detected as $hostname. Executing compilation for $hostname."

	# Source the environment setup file to enable module load
	. /vol0004/apps/oss/spack/share/spack/setup-env.sh

	# loading modules of hdf5 and fftw
	spack load /yhazdvl
	spack load /upvlzyl

	export GCC_FLAGS="-O3 -Kfast -Ksimd -funroll-loops -mcmodel=small -march=native -Kalign_commons -Kalign_loops -ffp-contract=fast"

	# Set CoToCoA mode by marco definition.
	export GCC_FLAGS="$GCC_FLAGS -D$cotocoa_mode"

	# Set mpi compiler.
	export CC="mpifccpx"
	export CFLAGS="$GCC_FLAGS"
	export FC="mpifrtpx"
	export FFLAGS="$GCC_FLAGS -Nalloc_assign -Dssurf=3 -lhdf5_fortran -lhdf5 -I /vol0004/apps/oss/spack-v0.21/opt/spack/linux-rhel8-a64fx/fj-4.10.0/hdf5-1.14.3-yhazdvld6vknkhmbcqrbl34ifsac2hao/include -I /vol0004/apps/oss/spack-v0.21/opt/spack/linux-rhel8-a64fx/fj-4.11.1/fftw-3.3.10-upvlzylw3inllggvpvxppuar7copfqdg/include"
	export FFLAGS="$FFLAGS -I./lib/cotocoa/include -L./lib/cotocoa/lib"

elif [[ $hostname == camphor* ]]; then # on camphor
    hostname=camphor
	echo "The current environment has been detected as $hostname. Executing compilation for $hostname."

	module load hdf5/1.12.2_intel-2022.3-impi
	module load fftw/3.3.10_intel-2022.3-impi

	export GCC_FLAGS="-ipo -O3 -no-prec-div -static_mpi -xHost -fp-model precise "

	# Set CoToCoA mode by marco definition.
	export GCC_FLAGS="$GCC_FLAGS -D$cotocoa_mode"

	# Set mpi compiler.
	export CC="mpiicc"
	export CFLAGS="$GCC_FLAGS -traceback -mcmodel=medium -shared-intel"

	export FC="mpiifort"
	export FFLAGS="$GCC_FLAGS -fpp -Dssurf=3 -traceback -mcmodel=medium -shared-intel"
	export FFLAGS="$FFLAGS -I./lib/cotocoa/include -L./lib/cotocoa/lib"
else
	echo This script is not compatible with the current environment. Exiting.
    exit 1
fi

# Install Fpm (Fortran package manager) if not exists.
PACKAGE="fpm"
if ! pip3 show "$PACKAGE" > /dev/null 2>&1; then
	echo "$PACKAGE is not currently installed. Initiating installation..."
	pip3 install --user "$PACKAGE"
else
	echo "$PACKAGE:ok"
fi

# Add path to use Fpm.
export PATH="$PATH:$HOME/.local/bin"

# Compile sub-code written by C lang.
make -C lib/mtarm
make -C lib/psort
make -C lib/ohhelp/ohhelp-1.1.1
make -C lib/cotocoa/src -f Makefile.env

export FPM_FC="$FC"
export FPM_CC="$CC"
export FPM_CXX="mpiixx"

export FPM_FFLAGS="$FFLAGS"
export FPM_FFLAGS="$FPM_FFLAGS -Llib/mtarm/ -Ilib/mtarm"
export FPM_FFLAGS="$FPM_FFLAGS -Llib/psort/ -Ilib/psort"
export FPM_FFLAGS="$FPM_FFLAGS -Llib/ohhelp/ohhelp-1.1.1/ -I./lib/ohhelp/ohhelp-1.1.1/"
#export FPM_FCFLAGS="$FPM_FCFLAGS -I$PWD/lib/cotocoa/include"
#export FPM_LDFLAGS="$FPM_LDFLAGS -L$PWD/lib/cotocoa/lib -lctca -lctca_f -lctca_mod"
export FPM_CFLAGS="$CFLAGS"

export CTCA_DEBUG=1
export CTCA_LOG_LEVEL=debug
export CTCA_LOG_FILE=ctca.log

export FPM_LDFLAGS=

fpm install --prefix ./
