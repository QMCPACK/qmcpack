#!/bin/bash

if [ -z $OMP_NUM_THREADS ]; then
  export OMP_NUM_THREADS=1
fi

if [ $# -lt 2 ] ; then
  echo "Need at least 2 arguments"
  exit 1
fi

ARGS=""

for i in `seq 3 $#`
do
  eval "arg=\${$i}"
  ARGS=$ARGS" "$arg
done

echo "COBALT_PARTNAME " $COBALT_PARTNAME
echo "MPI tasks       " $1
echo "Threads         " $OMP_NUM_THREADS
echo "Binary          " $2
echo "Input files     " $ARGS

runjob --np $1 -p 1 --block $COBALT_PARTNAME --envs OMP_NUM_THREADS=$OMP_NUM_THREADS : $2 $ARGS
