echo "COBALT_PARTNAME " $COBALT_PARTNAME
echo "MPI tasks       " $1
echo "Threads         " $OMP_NUM_THREADS
echo "Binary          " $2
echo "Input files     " $3 $4 $5 $6 $7

runjob --np $1 -p 1 --block $COBALT_PARTNAME --envs OMP_NUM_THREADS=$OMP_NUM_THREADS : $2 $3 $4 $5 $6 $7
