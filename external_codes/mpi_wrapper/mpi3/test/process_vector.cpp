#if COMPILATION_INSTRUCTIONS
mpicxx -O3 -std=c++14 `#-Wfatal-errors` $0 -o $0x.x -lboost_serialization -lboost_container && time mpirun -np 2 $0x.x $@ && rm -f $0x.x; exit
#endif

#include "alf/boost/mpi3/main.hpp"
//#include "alf/boost/mpi3/communicator.hpp"
#include "alf/boost/mpi3/process.hpp"
#include <boost/serialization/vector.hpp>
//#include "alf/boost/mpi3/detail/package_archive.hpp"

namespace boost{
namespace mpi3{

#if 0
template<class T, class A, typename = detail::basic_datatype<T> >
mpi3::process& operator<<(mpi3::process& p, std::vector<T, A> const& v){
	assert(v.size() <= std::size_t(std::numeric_limits<int>::max()));
	p.comm_.send_n(v.data(), v.size(), p.rank_);
	return p;
}

template<class T, class A, typename = detail::basic_datatype<T> >
mpi3::process& operator>>(mpi3::process&& p, std::vector<T, A>& v){
	v.resize(p.comm_.probe(p.rank_).count<T>());
	p.comm_.receive(v.begin(), v.end(), p.rank_);
	return p;
}

template<class T, class A, typename = detail::basic_datatype<T> > //class Vector, typename = decltype(Vector{}.data(), Vector{}.size())>
mpi3::communicator& operator>>(mpi3::communicator& comm, std::vector<T, A>& v){
	v.resize(comm.probe().count<T>());
	comm.receive(v.begin(), v.end());
}
#endif

}}

namespace mpi3 = boost::mpi3;
using std::cout;

struct long_long{
	long long value;
	long_long& operator=(long long v){value = v;}
};

template<class Archive>
void serialize(Archive& ar, long_long& l, unsigned = 0){
	ar & l.value;
}

int mpi3::main(int argc, char* argv[], mpi3::communicator& world){

	assert(world.size() == 2);
	long long size = 100000000;
	switch(world.rank()){
		case 0:{
			std::vector<long_long> v(size);
			std::iota(v.begin(), v.end(), 0.);
		//	assert(std::accumulate(v.begin(), v.end(), 0.) == size*(size-1)/2 );
			world[1] << v;
		}break;
		case 1:{
			std::vector<long_long> w;
			world[0] >> w;
			assert(w.size() == size);
			assert(w[45].value = 45);
		//	assert(std::accumulate(w.begin(), w.end(), 0.) == size*(size-1)/2 );
		}break;
	}
	return 0;

}

