#if COMPILATION_INSTRUCTIONS
mpicxx -O3 -std=c++17 `#-Wfatal-errors` $0 -o $0x.x && time mpirun -np 4 $0x.x $@ && rm -f $0x.x; exit
#endif

// based on https://computing.llnl.gov/tutorials/mpi/samples/C/array_mm.c
#include "alf/boost/mpi3/main.hpp"
#include "alf/boost/mpi3/version.hpp"
#include "alf/boost/mpi3/processor_name.hpp"
#include "alf/boost/multi/array.hpp"

#include<iostream>

using std::cout;
using std::endl;

int boost::mpi3::main(int argc, char* argv[], boost::mpi3::communicator const& world){
	assert( world.size() > 1 );

	if(world.master()){

		std::string name = "master";
		std::list<double> l = {1.,2.,3.,4.,5.,6.};

		boost::multi::array<double, 2> a({{62, 15}});
		boost::multi::array<double, 2> b({{15,7}});

		for(int i : a.extensions()[0])
			for(int j : a.extensions()[1])
				a[i][j] = i + j;

		for(int i : b.extensions()[0])
			for(int j : b.extensions()[1])
				b[i][j] = i*j;

		boost::multi::array<double, 2> c_serial({{62,7}});

		for(int i : c_serial.extensions()[0]){
			for(int j : c_serial.extensions()[1]){
				c_serial[i][j] = 0;
				for(int k : a.extensions()[1]){
					c_serial[i][j] += a[i][k] * b[k][j];
				}
			}
		}

		auto averow = a.extensions()[0].size()/(world.size()-1);
		auto extra  = a.extensions()[0].size()%(world.size()-1);
		
		auto offset = 0;
		for(int dest = 1; dest != world.size(); ++dest){
			auto p = world[dest];
			p << name;
			world.send_n(l.begin(), l.size(), dest);
			int rows = (dest < extra + 1)?(averow + 1): averow;
			p 
				<< offset 
				<< rows
			;
			world.send_n(&a[offset][0], rows*a.extensions()[1].size(), dest);
			world.send_n(&b[0][0], b.num_elements(), dest);
			offset += rows;
		}

		boost::multi::array<double, 2> c({{62,7}});
		for(int source = 1; source != world.size(); ++source){
			auto proc = world[source];
			int rows;
			proc 
				>> offset
				>> rows
			;
			world.receive_n(&c[offset][0], rows*c.extensions()[1].size(), source);
		}
		for(int i : c.extensions()[0])
			for(int j : c.extensions()[1])
				assert( c[i][j] - c_serial[i][j] == 0 );
	}else{
		auto proc = world[0];
		std::string name;
		proc >> name;
		std::list<double> l(6);
		world.receive_n(l.begin(), l.size(), 0);
		assert(( l == std::list<double>{1.,2.,3.,4.,5.,6.} ));
		int offset;
		int rows;
		proc >> offset >> rows;
		boost::multi::array<double, 2> a_partial({{rows, 15}});
		world.receive_n(a_partial.data(), a_partial.num_elements(), 0);
		boost::multi::array<double, 2> b({{15,7}});
		world.receive_n(b.data(), b.num_elements(), 0);

		boost::multi::array<double, 2> c_partial({{rows,7}});
		for(int i = 0; i != c_partial.extensions()[0].size(); ++i){
			for(int j = 0; j != c_partial.extensions()[1].size(); ++j){
				c_partial[i][j] = 0;
				for(int k = 0; k != a_partial.extensions()[1].size(); ++k){
					c_partial[i][j] += a_partial[i][k] * b[k][j];
				}
			}
		}
		proc 
			>> offset 
			>> rows
		;
		world.send_n(c_partial.data(), c_partial.num_elements(), 0);
	}
	cout << "ok" << endl;
}

