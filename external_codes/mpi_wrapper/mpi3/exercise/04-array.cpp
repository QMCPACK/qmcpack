#if COMPILATION_INSTRUCTIONS
mpicxx -O3 -std=c++17 `#-Wfatal-errors` $0 -o $0x.x && time mpirun -np 4 $0x.x $@ && rm -f $0x.x; exit
#endif
// based on https://computing.llnl.gov/tutorials/mpi/samples/C/mpi_array.c
#include "alf/boost/mpi3/main.hpp"
#include "alf/boost/mpi3/version.hpp"
#include "alf/boost/mpi3/processor_name.hpp"

#include<iostream>

using std::cout;
using std::endl;

double update(std::vector<double>& data, int offset, int chunk_size){
	double sum = 0;
	for(int i = offset; i != offset + chunk_size; ++i){
		data[i] += i*1.0;
		sum += data[i];
	}
	return sum;
}

int boost::mpi3::main(int argc, char* argv[], boost::mpi3::communicator const& world){
	if(world.size() % 4 != 0){
		if(world.master()) cout << "Quitting. Need an even number of tasks: numtasks = " << world.size() << std::endl;
		return 1;
	}

	int const array_size = 16000000;
	int chunk_size = array_size/world.size();

	if(world.master()){
		std::vector<double> data(array_size);
		{
			double sum = 0;
			for(int i = 0; i != data.size(); ++i){
				data[i] = i*1.0;
				sum += data[i];
			}
			cout << "serial sum " << sum << endl;
		}

		int offset = chunk_size;
		for(int dest = 1; dest != world.size(); ++dest){
			world.send_n(data.begin() + offset, chunk_size, dest);
			offset += chunk_size;
		}

		double partial_sum = 0;
		for(int i = 0; i != chunk_size; ++i)
			partial_sum += data[i];

		double parallel_sum = world.reduce(partial_sum);
		cout << "parallel sum is " << parallel_sum << endl;
	}else{
		std::vector<double> partial_data(chunk_size);
		world.receive_n(partial_data.begin(), chunk_size, 0);

		double partial_sum = 0;
		for(int i = 0; i != partial_data.size(); ++i)
			partial_sum += partial_data[i];

		world.reduce(partial_sum);
	}


}

