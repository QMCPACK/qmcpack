#if COMPILATION_INSTRUCTIONS
mpic++ -O3 -std=c++14 -Wall -Wextra -D_MAKE_BOOST_SERIALIZATION_HEADER_ONLY `#-lboost_serialization` $0 -o $0x.x && time mpirun -n 3 $0x.x $@ && rm -f $0x.x; exit
#endif

#include "../../mpi3/main.hpp"
#include "../../mpi3/communicator.hpp"
#include "../../mpi3/process.hpp"


namespace mpi3 = boost::mpi3;
using std::cout;

// nontrivial nonpod class
struct B{
	std::string name_ = "unnamed"; 
	int n_ = 0;
	double* data = nullptr;
	B() = default;
	B(int n) : n_(n), data(new double[n]){
		for(int i = 0; i != n_; ++i) data[i] = 0.;
	}
	B(B const& other) : name_(other.name_), n_(other.n_), data(new double[other.n_]){}
	B& operator=(B const& other){
		name_ = other.name_;
		n_ = other.n_; 
		delete[] data; 
		data = new double[other.n_];
		for(int i = 0; i != n_; ++i) data[i] = other.data[i];
		return *this;
	}
	~B(){delete[] data;}
};

// nonintrusive serialization
template<class Archive>
void save(Archive & ar, B const& b, const unsigned int){
	ar << b.name_ << b.n_ << boost::serialization::make_array(b.data, b.n_);
}
template<class Archive>
void load(Archive & ar, B& b, const unsigned int){
	ar >> b.name_ >> b.n_;
	delete[] b.data; b.data = new double[b.n_];
	ar >> boost::serialization::make_array(b.data, b.n_);
}
BOOST_SERIALIZATION_SPLIT_FREE(B)

int mpi3::main(int, char*[], mpi3::communicator world){
	assert( world.size() > 1 );

	switch(world.rank()){
		case 0 : {
			std::vector<B> v(5, B(3));
			v[2].data[2] = 3.14;
			world.send(v.begin(), v.end(), 1, 123);
		}; break;
		case 1 : {
			std::vector<B> v(5);
			world.receive(v.begin(), v.end(), 0, 123);
			assert(v[2].data[2] == 3.14);
		}; break;
	}
	switch(world.rank()){
		case 0 : {
			B b1(4); b1.data[2] = 4.5;
			world[1] << b1;
		}; break;
		case 1 : {
			B b2;
			world[0] >> b2;
			assert( b2.data[2] == 4.5 );
		}; break;
	}
	
	return 0;
}

