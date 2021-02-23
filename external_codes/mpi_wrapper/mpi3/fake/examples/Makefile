all : abort.x MPI_Send.x MPI_Recv.x MPI_Comm_dup.x
	mpirun -np 1 ./abort.x
	mpirun -np 1 ./MPI_Send.x
	mpirun -np 1 ./MPI_Recv.x
	mpirun -np 1 ./MPI_Comm_dup.x
	rm -f $^

abort.x : abort.cpp
	mpic++ abort.cpp -o abort.x

MPI_Send.x : MPI_Send.cpp
	mpic++ MPI_Send.cpp -o MPI_Send.x

MPI_Recv.x : MPI_Recv.cpp
	mpic++ MPI_Recv.cpp -o MPI_Recv.x

MPI_Comm_dup.x : MPI_Comm_dup.cpp
	mpic++ MPI_Comm_dup.cpp -o MPI_Comm_dup.x

clean:
	rm -f *.x

load:
	module unload mpi
	
