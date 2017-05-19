#!/bin/bash

module load env/dev

NMPI=2
NTHREAD=8

export OMP_NUM_THREADS=$NTHREAD
BIN=/home/yyang173/soft/qmcpack/build/bin/qmcpack
mpirun -np $NMPI $BIN vmc.xml >& vmc.out

./latdev_check.py
