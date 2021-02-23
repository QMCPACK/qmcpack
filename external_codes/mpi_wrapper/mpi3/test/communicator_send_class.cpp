#if COMPILATION_INSTRUCTIONS
mpic++ -O3 -std=c++14 -Wfatal-errors -D_MAKE_BOOST_SERIALIZATION_HEADER_ONLY `#-lboost_serialization` $0 -o $0x.x && time mpirun -n 2 $0x.x $@ && rm -f $0x.x; exit
#endif
//  (C) Copyright Alfredo A. Correa 2018.

#include "../../mpi3/main.hpp"
#include "../../mpi3/communicator.hpp"
//#include "../../mpi3/detail/package_archive.hpp"

#include<boost/serialization/vector.hpp>
#include<boost/serialization/utility.hpp> // serialize std::pair
#include<set>

namespace mpi3 = boost::mpi3;
using std::cout;

struct A{
	std::string name_ = "unnamed"; 
	int n_ = 0;
	double* data = nullptr;
	A() = default;
	A(int n) : n_(n), data(new double[n]){}
	A(A const& other) : name_(other.name_), n_(other.n_), data(new double[other.n_]){}
	A& operator=(A const& other){
		name_ = other.name_;
		n_ = other.n_; 
		delete[] data; 
		data = new double[other.n_];
		for(int i = 0; i != n_; ++i) data[i] = other.data[i];
		return *this;
	}
	~A(){delete[] data;}
	// intrusive serialization
    template<class Archive>
    void save(Archive & ar, const unsigned int) const{
		ar << name_ << n_ << boost::serialization::make_array(data, n_);
    }
    template<class Archive>
    void load(Archive & ar, const unsigned int){
		ar >> name_ >> n_;
		delete[] data; data = new double[n_];
		ar >> boost::serialization::make_array(data, n_);
    }
    BOOST_SERIALIZATION_SPLIT_MEMBER()
};

struct B{
	std::string name_ = "unnamed"; 
	int n_ = 0;
	double* data = nullptr;
	B() = default;
	B(int n) : n_(n), data(new double[n]){}
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
			std::vector<std::vector<double>> buffer(10, std::vector<double>(20));
			buffer[4][5] = 6.1;
			world.send(buffer.begin(), buffer.end(), 1, 123);
		}; break;
		case 1 : {
			std::vector<std::vector<double>> buffer(10, std::vector<double>(20));
			world.receive(buffer.begin(), buffer.end(), 0, 123);
			assert( buffer[4][5] == 6.1 );
		}; break;
	}
	switch(world.rank()){
		case 0 : {
			std::vector<double> buffer(10);
			iota(begin(buffer), end(buffer), 0);
			world.send(begin(buffer), end(buffer), 1, 123);
		}; break;
		case 1 : {
			std::vector<double> v(10);
		//	world.receive(std::back_inserter(v), 0, 123);
			auto it = world.receive(begin(v), 0, 123);
			assert(it == v.end() and v[3] == 3.);
		}; break;
	}
	switch(world.rank()){
		case 0: {
			std::map<int, std::vector<double>> m;
			m[2] = std::vector<double>(2);
			m[5] = std::vector<double>(5);
			world.send(begin(m), end(m), 1, 123);
		}; break;
		case 1: {
			std::vector<std::pair<int, std::vector<double>>> v(2);
			world.receive(begin(v), end(v), 0, 123);
			assert(( v[1] == std::pair<int, std::vector<double>>{5, std::vector<double>(5)} ));
		}; break;
	}
	switch(world.rank()){
		case 0 : {
			std::vector<A> v(5, A(3));
			v[2].data[2] = 3.14;
			world.send(begin(v), end(v), 1, 123);
		}; break;
		case 1 : {
			std::vector<A> v(5);
			world.receive(begin(v), end(v), 0, 123);
			assert(v[2].data[2] == 3.14);
		}; break;
	}
	switch(world.rank()){
		case 0 : {
			std::vector<B> v(5, B(3));
			v[2].data[2] = 3.14;
			world.send(begin(v), end(v), 1, 123);
		}; break;
		case 1 : {
			std::vector<B> v(5);
			world.receive(begin(v), end(v), 0, 123);
			assert(v[2].data[2] == 3.14);
		}; break;
	}

	return 0;
}

