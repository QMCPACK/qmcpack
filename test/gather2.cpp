#include "../../mpi3/main.hpp"
#include "../../mpi3/communicator.hpp"
#include "../../mpi3/detail/iterator.hpp"

namespace mpi3 = boost::mpi3;

auto mpi3::main(int/*argc*/, char**/*argv*/, mpi3::communicator world) -> int try{

	using dd = std::tuple<double, double>;

	std::vector<dd> v_local(10, dd{world.rank(), world.rank() + 1});
	std::vector<dd> v(world.root()?v_local.size()*static_cast<std::size_t>(world.size()):0);
	auto last = world.gather(begin(v_local), end(v_local), begin(v));
	if(world.root()){
		assert(last == end(v)); 
		assert(( v[ 0] == dd{0.,1.} ));
		assert(( v[10] == dd{1.,2.} ));
		assert(( v[20] == dd{2.,3.} ));
	}else{
		assert( last == begin(v) );
	}
	return 0;
} catch(...) {return 1;}
