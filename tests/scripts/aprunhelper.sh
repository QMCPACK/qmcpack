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

echo "MPI tasks       " $1
echo "Threads         " $OMP_NUM_THREADS
echo "Binary          " $2
echo "Input files     " $ARGS

aprun -n $1 -N $1 -d $OMP_NUM_THREADS -cc depth -j 1 $2 $ARGS
