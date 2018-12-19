#if COMPILATION_INSTRUCTIONS
time mpicxx -O3 -std=c++14 -Wfatal-errors -Wall $0 -o $0x.x -lboost_system && time mpirun -np 1 $0x.x $@ && rm -f $0x.x; exit
#endif

#include "alf/boost/mpi3/shm/managed_shared_memory.hpp" // there is a bug in boost 1.64, this needs to be included first
#include "alf/boost/mpi3/main.hpp"

#include<atomic>
#include<boost/interprocess/offset_ptr.hpp>

namespace mpi3 = boost::mpi3;
using std::cout; 

struct list_node{
//	boost::interprocess::offset_ptr<list_node> next;
	list_node* next;
	int value;
};

int mpi3::main(int argc, char* argv[], mpi3::communicator& world){

	mpi3::communicator node = world.split_shared(0);

	mpi3::shm::managed_shared_memory segment(node, 65536);
	list_node* prev = 0;
	list_node* current;
	list_node* first;

	list_node* ln = static_cast<list_node*>(segment.allocate(sizeof(list_node)));
	assert(ln != nullptr);
	ln -> value = 0;
	ln -> next = 0;
	list_node* ln2 = static_cast<list_node*>(segment.allocate(sizeof(list_node)));
	assert(ln2 != nullptr);
	ln2 -> value = 1;
//	ln2 -> next = ln;

	for(int i = 0; i < 10; ++i, prev = current){
		current = new list_node;
	//	current = static_cast<list_node*>(segment.allocate(sizeof(list_node)));
		current->value = i;
		current->next  = 0;
		if(!prev) first = current;
		else prev->next = current;
	}

#if 0
	for(current = first; current; ){
		prev = current;
		current = current->next;
		segment.deallocate(prev);//.get());
	}
#endif

	return 0;
}

