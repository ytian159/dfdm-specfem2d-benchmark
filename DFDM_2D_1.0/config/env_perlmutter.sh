#!/bin/bash

module load cmake
module load PrgEnv-cray
module load gcc
export OPEN_BLAS_CMAKE_PATH=/pscratch/sd/m/mgawan/openblas_install/OpenBLAS/install/lib/cmake


# export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=16
# export OMP_PROC_BIND=close
# export OMP_PLACES=cores