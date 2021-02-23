#if COMPILATION_INSTRUCTIONS
mpic++ $0 -o $0x -lboost_serialization&&mpirun -n 3 $0x&&rm $0x;exit
#endif

#include "../../mpi3/main.hpp"
#include "../../mpi3/communicator.hpp"
#include "../../mpi3/detail/iterator.hpp"
#include<numeric> //iota

namespace mpi3 = boost::mpi3;
using std::cout;

template<class T> void f(int);

int mpi3::main(int, char*[], mpi3::communicator world){
	assert( world.size() > 2 );
{
	using T = std::tuple<double, double>;
	std::vector<T> v_local(10, T{world.rank(), world.rank()});
	std::vector<T> v(v_local.size()*world.size());
	auto end = world.all_gather_n(v_local.begin(), v_local.size(), v.begin());
	assert(end == v.end());
	assert(( v[ 0] == T{0.,0.} ));
	assert(( v[ 9] == T{0.,0.} ));
	assert(( v[10] == T{1.,1.} ));
	assert(( v[19] == T{1.,1.} ));
	assert(( v[20] == T{2.,2.} ));
	assert(( v[29] == T{2.,2.} ));
}
{
	using T = std::tuple<double, double>;
	std::vector<T> v_local(world.rank() + 5, T{world.rank(), world.rank()});
	std::vector<T> v(1000, T{-99., -99.});
	auto d_last = world.all_gatherv_n(begin(v_local), v_local.size(), begin(v));
	
	int predict_size = 0;
	for(auto i = 0; i != world.size(); ++i) predict_size += i + 5;
	assert( std::distance(begin(v), d_last) == predict_size );
	
	if(world.rank()==1){
		cout<< std::distance(begin(v), d_last) <<std::endl;
		for(auto it = begin(v); it != d_last; ++it)
			cout<<"("<<std::get<0>(*it)<<' '<<std::get<1>(*it)<<"), "<<std::endl;
	}
	assert(( v[ 0] == T{0.,0.} ));
	assert(( v[ 4] == T{0.,0.} ));
	assert(( v[ 5] == T{1.,1.} ));
	assert(( v[10] == T{1.,1.} ));
	assert(( v[11] == T{2.,2.} ));
	assert(( v[17] == T{2.,2.} ));
}
{
	using T = std::tuple<double, double>;
	std::vector<T> v_local(world.rank() + 5, T{world.rank(), world.rank()});
	std::vector<T> v(1000, T{-99., -99.});
	auto d_last = world.all_gatherv(begin(v_local), end(v_local), begin(v));
	assert( d_last < end(v) );
}
{
// initialize data
	using T = double;
	assert( world.size() == 3 );
	std::vector<T> v_loc;
	switch(world.rank()){
		case 0: v_loc = {0., 0., 0.}        ; break;
		case 1: v_loc = {1., 1., 1., 1.}    ; break;
		case 2: v_loc = {2., 2., 2., 2., 2.}; break;
	}
// gather communication
	std::vector<T> v;
//	v.reserve(v_local.size()*world.size()); // to avoid _some_ reallocations
	world.all_gatherv(begin(v_loc), end(v_loc), std::back_inserter(v)); 
//	v.shrink_to_fit(); // to save _some_ memory

// check communication
	assert((v==std::vector<T>{0., 0., 0., 1., 1., 1., 1., 2., 2., 2., 2., 2.}));
}
	return 0;
}

