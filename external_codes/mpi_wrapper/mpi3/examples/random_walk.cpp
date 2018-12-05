#if COMPILATION_INSTRUCTIONS
mpicxx -O3 -std=c++14 -Wall `#-Wfatal-errors` $0 -o $0x.x -lboost_serialization && time mpirun -np 4 $0x.x $@ && rm -f $0x.x; exit
#endif

#include "alf/boost/mpi3/main.hpp"
#include "alf/boost/mpi3/communicator.hpp"
#include "alf/boost/mpi3/process.hpp"
#include "alf/boost/mpi3/detail/package_archive.hpp"

#include<boost/serialization/vector.hpp>

#include <random>

namespace mpi3 = boost::mpi3;
using std::cout;

struct walker{
	int location;
	int steps_left;
	template<class Archive>
	void serialize(Archive& ar, const unsigned int){
		ar & location & steps_left;
	}
};

int mpi3::main(int, char*[], mpi3::communicator& world){
	int domain_size = 20;
	int max_walk_size = 100;
	int num_walkers_per_proc = 10;

	// decompose domain
	assert(world.size() <= domain_size);
	int subdomain_start = domain_size/world.size()*world.rank();
	int subdomain_size  = domain_size/world.size();
	if(world.rank() == world.size() - 1) subdomain_size += domain_size%world.size();

	// initialize walkers
	std::vector<walker> incomming_walkers(num_walkers_per_proc);
	std::vector<walker> outgoing_walkers;

	std::random_device rd;
	std::mt19937 gen(rd()*world.rank());
    std::uniform_int_distribution<> dis(0, max_walk_size);
	std::generate(
		incomming_walkers.begin(), incomming_walkers.end(), 
		[&]{return walker{subdomain_start, dis(gen)};}
	);

	cout << "process " << world.rank() << " initialized " << num_walkers_per_proc << " walkers in subdomain " << subdomain_start << " - " << subdomain_start + subdomain_size - 1 << std::endl;
	
	// determine the maximum amount of sends and receives needed to complete all walkers
	int maximum_sends_recvs = max_walk_size/(domain_size/world.size()) + 1;
	for(int m = 0; m != maximum_sends_recvs; ++m){
		// process all incomming_walkers
		for(auto& w : incomming_walkers){
			while(w.steps_left){
				if(w.location == subdomain_start + subdomain_size){
					if(w.location == domain_size) w.location = 0;
					outgoing_walkers.push_back(w);
					break;
				}else{
					--w.steps_left;
					++w.location;
				}
			}
		}
		cout << "process " << world.rank() << " sending " << outgoing_walkers.size() << " outgoing walkers to process " << (world.rank() + 1)%world.size() << std::endl;
		if(world.rank()%2 == 0){
			world[(world.rank() + 1)%world.size()] << outgoing_walkers;
			outgoing_walkers.clear();
			world[world.rank()?world.rank()-1:world.size()-1] >> incomming_walkers;
		}else{
			world[world.rank()?world.rank()-1:world.size()-1] >> incomming_walkers;
			world[(world.rank() + 1)%world.size()] << outgoing_walkers;
			outgoing_walkers.clear();
		}
		cout << "process " << world.rank() << " received " << incomming_walkers.size() << " incomming_walkers" << std::endl;
	}
	cout << "process " << world.rank() << " done" << std::endl;
	return 0;
}

