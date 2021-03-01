#if COMPILATION_INSTRUCTIONS
mpic++ $0 -o $0x&&mpirun -n 4 $0x $@&&rm $0x;exit
#endif

#include "../../mpi3/main.hpp"
#include "../../mpi3/shared_communicator.hpp"
#include "../../mpi3/ostream.hpp"

#include<iostream>

namespace mpi3 = boost::mpi3;
using std::cout;

void check_isa_communicator(mpi3::communicator const&){}
void check_isa_shared_comm (mpi3::shared_communicator const&){}

int mpi3::main(int, char*[], mpi3::communicator world){

	mpi3::ostream wout(world);

	auto node = world.split_shared();
	wout << "I am rank " << node.rank() << " in comm " << node.name() << std::endl;
	wout << "----" << std::endl;

	auto core = world.split_shared(mpi3::communicator_type::core);
	wout << "I am rank " << core.rank() << " in comm " << core.name() << std::endl;
	wout << "----" << std::endl;

	auto hw = world.split_shared(mpi3::communicator_type::hw_thread);
	wout << "I am rank " << hw.rank() << " in comm " << hw.name() << std::endl;
	wout << "----" << std::endl;

	auto l1 = world.split_shared(mpi3::communicator_type::l1_cache);
	wout << "I am rank " << l1.rank() << " in comm " << l1.name() << std::endl;
	wout << "----" << std::endl;

	auto l2 = world.split_shared(mpi3::communicator_type::l2_cache);
	wout << "I am rank " << l2.rank() << " in comm " << l2.name() << std::endl;
	wout << "----" << std::endl;

	auto l3 = world.split_shared(mpi3::communicator_type::l3_cache);
	wout << "I am rank " << l3.rank() << " in comm " << l3.name() << std::endl;
	wout << "----" << std::endl;

	auto socket = world.split_shared(mpi3::communicator_type::socket);
	wout << "I am rank " << socket.rank() << " in comm " << socket.name() << std::endl;
	wout << "----" << std::endl;

	auto numa = world.split_shared(mpi3::communicator_type::numa);
	wout << "I am rank " << numa.rank() << " in comm " << numa.name() << std::endl;
	wout << "----" << std::endl;

	auto board = world.split_shared(mpi3::communicator_type::board);
	wout << "I am rank " << board.rank() << " in comm " << board.name() << std::endl;
	wout << "----" << std::endl;

	auto host = world.split_shared(mpi3::communicator_type::host);
	wout << "I am rank " << host.rank() << " in comm " << host.name() << std::endl;
	wout << "----" << std::endl;

	auto cu = world.split_shared(mpi3::communicator_type::cu);
	wout << "I am rank " << cu.rank() << " in comm " << cu.name() << std::endl;
	wout << "----" << std::endl;

	return 0;
#if 0
	mpi3::shared_communicator node2 = node.split(node.rank() % 2, node.rank());

	wout << "I am rank " << node.rank() << " in comm " << node.name() << " and rank " << node2.rank() << " in " << node2.name() << std::endl;

	check_isa_shared_comm(node);
	check_isa_communicator(node);

	mpi3::communicator virtual_node{node};
	assert( &virtual_node != &node );

	mpi3::communicator virtual_node2;
	virtual_node2 = node2;
	assert( &virtual_node2 != &node2 );

	mpi3::communicator virtual_node3  = world.split_shared();
	assert( &virtual_node3 != &world );

	assert( sizeof(mpi3::communicator) == sizeof(mpi3::shared_communicator) );
#endif
	return 0;
}

