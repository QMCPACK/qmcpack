#if COMPILATION_INSTRUCTIONS
mpicxx.mpich -g $0 -o $0x -lboost_serialization&&mpirun.mpich -n 3 valgrind --error-exitcode=1345 $0x&&rm $0x;exit
#endif
// © Alfredo A. Correa 2018-2021

#include "../../mpi3/main.hpp"
#include "../../mpi3/communicator.hpp"
#include "../../mpi3/ostream.hpp"
#include "../../mpi3/process.hpp"

#include<boost/serialization/utility.hpp>

#include<list>

namespace mpi3 = boost::mpi3;

using std::cout;
using std::list;
using std::vector;

void part1(mpi3::communicator& world) {
	vector<std::pair<double, int>> small(10, {0.0, world.rank()});
	vector<std::pair<double, int>> large(world.root() ? small.size() * static_cast<std::size_t>(world.size()) : 0, std::pair<double, int>(0.0, -10));

	auto it = world.gather_n(small.begin(), small.size(), large.begin(), 0);
	assert(it == large.end());

	if(world.rank() == 0) {
		assert(it != large.begin());
		assert((large[9] == std::pair<double, int>(0.0, 0)));
		assert((large[11] == std::pair<double, int>(0.0, 1)));
	} else {
		assert(it == large.begin());
	}
}

void part2(mpi3::communicator& world) {
	vector<std::pair<double, int>> small(10, {0., world.rank()});
	vector<std::pair<double, int>> large(world.root() ? small.size() * static_cast<std::size_t>(world.size()) : 0, std::pair<double, int>(0., -1));

	auto it = world.gather_n(small.begin(), small.size(), large.begin());
	assert(it == large.end());

	if(world.root()) {
		assert(it != large.begin());
		assert((large[9] == std::pair<double, int>(0., 0)));
		assert((large[11] == std::pair<double, int>(0., 1)));
	} else {
		assert(it == large.begin());
	}
}

void part3(mpi3::communicator& world) {
	list<double> small(10, world.rank());
	vector<double> large(world.root()?small.size()*static_cast<std::size_t>(world.size()):0, -1.0);

	world.gather(small.begin(), small.end(), large.begin(), 0);
	if(world.root()) {
		cout << "large: ";
		for(auto& e : large) {cout<< e <<" ";}  // NOLINT(altera-unroll-loops)
		cout << '\n';
	}
	if(world.root()) {
		assert(large[ 1] == 0);
		assert(large[11] == 1);
		assert(large[21] == 2);
	}
}

void part4(mpi3::communicator& world) {
	auto val = std::string{"5.1 opa"};

	using T = decltype(val);

	vector<T> small(10, val);
	vector<T> large(world.root() ? small.size() * static_cast<std::size_t>(world.size()) : 0);

	world.gather(small.begin(), small.end(), large.begin(), 0);

	if(world.rank() == 0) {
		assert(all_of(large.begin(), large.end(), [val](auto& e) { return val == e; }));
	}
}

void part5(mpi3::communicator& world) {
	auto Lval  = std::to_string(world.rank() + 1000);
	auto vals0 = (world[0] |= Lval);
	if(world.rank() == 0) {
		assert(vals0.size() - static_cast<std::size_t>(world.size()) == 0);
		assert(vals0[2] == "1002");
	} else {
		assert(vals0.empty());
	}
}

auto mpi3::main(int/*argc*/, char**/*argv*/, mpi3::communicator world) -> int try{

	mpi3::vector<double> v;

	assert( world.size() > 2);

	part1(world);
	part2(world);
	part3(world);
	part4(world);
	part5(world);

	return 0;

}catch(...){
	return 1;
}
