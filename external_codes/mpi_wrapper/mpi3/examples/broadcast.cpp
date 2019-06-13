#if COMPILATION_INSTRUCTIONS
mpicxx -O3 -std=c++14 -Wall -Wfatal-errors $0 -o $0x.x -lboost_serialization && time mpirun -np 4 $0x.x $@ && rm -f $0x.x; exit
#endif

#include "alf/boost/mpi3/main.hpp"
#include "alf/boost/mpi3/process.hpp"
#include "alf/boost/mpi3/detail/package_archive.hpp"

#include<boost/serialization/vector.hpp>

namespace mpi3 = boost::mpi3;
using std::cout;

void my_broadcast(mpi3::communicator& comm, std::vector<double>& data, int root){
	if(comm.rank() == root){
		for(int i = 0; i != comm.size(); ++i){
			if(i != comm.rank()) comm[i] << data;
		}
	}else{
		comm[root] >> data;
	}
}

void my_bcast(void* data, int count, MPI_Datatype datatype, int root,
              MPI_Comm communicator) {
  int world_rank;
  MPI_Comm_rank(communicator, &world_rank);
  int world_size;
  MPI_Comm_size(communicator, &world_size);

  if (world_rank == root) {
    // If we are the root process, send our data to everyone
    int i;
    for (i = 0; i < world_size; i++) {
      if (i != world_rank) {
        MPI_Send(data, count, datatype, i, 0, communicator);
      }
    }
  } else {
    // If we are a receiver process, receive the data from the root
    MPI_Recv(data, count, datatype, root, 0, communicator, MPI_STATUS_IGNORE);
  }
}

int mpi3::main(int argc, char* argv[], mpi3::communicator& world){
	{
		std::vector<double> data;
		if(world.rank() == 0){data.resize(10); std::iota(data.begin(), data.end(), 0);}
		my_broadcast(world, data, 0);
		assert( data[5] == 5 );
	}
	// ... less efficient but same effect as ...
	{
		std::vector<double> data;
		if(world.rank() == 0){data.resize(10); std::iota(data.begin(), data.end(), 0);}
		world[0] & data;
		assert( data[5] == 5 );
	}
	int num_elements = 100000;
	int num_trials = 1000;
	{
		int world_rank;
		MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

		double total_my_bcast_time = 0.0;
		double total_mpi_bcast_time = 0.0;
		int i;
		double* data = (double*)malloc(sizeof(double) * num_elements);
		std::iota(data, data + num_elements, 0);
		assert(data != NULL);

		for (i = 0; i < num_trials; i++) {
			// Time my_bcast
			// Synchronize before starting timing
			MPI_Barrier(MPI_COMM_WORLD);
			total_my_bcast_time -= MPI_Wtime();
			my_bcast(data, num_elements, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			// Synchronize again before obtaining final time
			MPI_Barrier(MPI_COMM_WORLD);
			total_my_bcast_time += MPI_Wtime();

			// Time MPI_Bcast
			MPI_Barrier(MPI_COMM_WORLD);
			total_mpi_bcast_time -= MPI_Wtime();
			MPI_Bcast(data, num_elements, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			MPI_Barrier(MPI_COMM_WORLD);
			total_mpi_bcast_time += MPI_Wtime();
		}

		// Print off timing information
		if (world_rank == 0) {
			printf("Data size = %d, Trials = %d\n", num_elements * (int)sizeof(double),
				   num_trials);
			printf("Avg my_bcast time = %lf\n", total_my_bcast_time / num_trials);
			printf("Avg MPI_Bcast time = %lf\n", total_mpi_bcast_time / num_trials);
		}
		free(data);
	}
	{

		std::vector<double> data(num_elements); std::iota(data.begin(), data.end(), 0);

		double total_my_broadcast_time = 0.0;
		double total_broadcast_time = 0.0;
		for(int i = 0; i != num_trials; ++i){
		//	if(world.rank() == 0) cout << i << std::endl;
			world.barrier();
			total_my_broadcast_time -= mpi3::wall_time();
			my_broadcast(world, data, 0);
			world.barrier();
			total_my_broadcast_time += mpi3::wall_time();
			world.barrier();
			total_broadcast_time -= mpi3::wall_time();
		//	world.broadcast(data.begin(), data.end(), 0);
		//	MPI_Bcast(data.data(), num_elements, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		//	world[0] & data;
			int n = data.size(); world[0] & n; data.resize(n);
			world.broadcast(data.begin(), data.end(), 0);
			world.barrier();
			total_broadcast_time += mpi3::wall_time();
		}
		if(world.rank() == 0){
			cout << "data size = " << num_elements * (int)sizeof(double) << ", trials = " << num_trials << '\n';
			cout << "avg my_broadcast time = " << total_my_broadcast_time / num_trials << '\n';
			cout << "avg broadcast time = " << total_broadcast_time / num_trials << '\n';
		}
	}

	return 0;
}

