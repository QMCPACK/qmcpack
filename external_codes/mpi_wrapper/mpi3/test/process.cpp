#if COMPILATION_INSTRUCTIONS
mpicxx -O3 -std=c++14 `#-Wfatal-errors` $0 -o $0x.x -lboost_serialization && time mpirun -np 4 $0x.x $@ && rm -f $0x.x; exit
#endif

#include "../../mpi3/main.hpp"
#include "../../mpi3/communicator.hpp"
#include "../../mpi3/process.hpp"

//#include "../../mpi3/detail/package_archive.hpp"

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

namespace mpi3 = boost::mpi3;
using std::cout;

int mpi3::main(int, char*[], mpi3::communicator world){

	assert(world.size() > 1);

	if(world.rank() == 0){
		int a = 5;
		world[1] << a;
	}else if(world.rank() == 1){
		int a = -1;
		world[0] >> a; // specific source (any tag)
		assert(a == 5);
	}

	if(world.rank() == 0){
		int a = 7;
		world[1] << a;
	}else if(world.rank() == 1){
		int a = -1;
		world >> a; // any source (any tag)
		assert(a == 7);
	}

	int b = world.rank();
	world[2] & b; // broadcast (from rank 2)
	assert( b == 2 );

	return 0;

	if(world.root()){
		B b1(4); b1.data[2] = 4.5;
		world[1] << b1;
	}else{
		B b2;
		world[0] >> b2;
		assert( b2.data[2] == 4.5 );
	}
	return 0;
}

