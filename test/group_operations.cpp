#if COMPILATION_INSTRUCTIONS
mpicxx -O3 -std=c++14 `#-Wfatal-errors` $0 -o $0x.x && time mpirun -np 8s $0x.x $@ && rm -f $0x.x; exit
#endif

#include "../../mpi3/environment.hpp"
#include "../../mpi3/group.hpp"
#include "../../mpi3/communicator.hpp"

using std::cout;
namespace mpi3 = boost::mpi3;

int main(int argc, char* argv[]){
	mpi3::environment env(argc, argv);

	mpi3::communicator world = env.world();
	mpi3::communicator comm = env.world();

	mpi3::group basegroup(world);
	int grp_rank = basegroup.rank();
	assert(basegroup.rank() == comm.rank());

    /* Form a new communicator with inverted ranking */
	mpi3::communicator newcomm = comm.split(0, comm.size() - comm.rank());
	mpi3::group g1(newcomm);
	std::vector<int> ranks_out = g1.translate_ranks(basegroup);
	for(int i = 0; i != ranks_out.size(); ++i) assert(ranks_out[i] == comm.size() - 1 - i);

	assert( basegroup.compare(g1) == mpi3::similar );
	mpi3::communicator dupcomm = comm;
	mpi3::group g2(dupcomm);
	assert( basegroup.compare(g2) == mpi3::identical );

	mpi3::communicator splitcomm = comm.split(comm.rank() < comm.size()/2, comm.rank());
	mpi3::group g3(splitcomm);
	assert( compare(basegroup, splitcomm) == mpi3::unequal );
	
	mpi3::group g3a = basegroup.include({comm.rank(), (comm.rank() + 1)%comm.size()});
	mpi3::group g3b = basegroup.include({comm.rank(), (comm.rank() + comm.size() -1)%comm.size()});
	assert( compare(g3a, g3b) == mpi3::unequal );


	std::vector<int> ranks(basegroup.size()); 
	std::iota(ranks.begin(), ranks.end(), 0);
	mpi3::group g4 = basegroup.exclude_n(ranks.begin(), 1);
	mpi3::group g5 = basegroup.exclude_n(ranks.begin() + 1, ranks.size() - 1);
	mpi3::group g6 = set_union(g5, g4); // order matters
	assert( basegroup == g6 );

	mpi3::group g7 = set_union(basegroup, g4);
	assert( basegroup == g7 );

	mpi3::group g8 = basegroup.range_exclude({{1, basegroup.size()-1, 1}});
	assert( g5 == g8 );

	mpi3::group g9 = set_intersection(basegroup, g4);
	assert( g9 == g4 );

	mpi3::group g10 = basegroup.range_exclude({{0, basegroup.size()-1, 1}});
//	assert( g10 == mpi3::group::empty() ); // g10 == mpi3::group()
//	assert( g10.size() == 0 );
	assert( g10.empty() );
}

