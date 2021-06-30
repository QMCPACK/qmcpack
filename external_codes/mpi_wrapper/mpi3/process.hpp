#ifdef COMPILATION_INSTRUCTIONS /* -*- indent-tabs-mode: t -*- */
(echo "#include\""$0"\"" > $0x.cpp) && mpic++ -O3 -std=c++14 -Wall -D_TEST_BOOST_MPI3_PROCESS $0x.cpp -o $0x.x -lboost_serialization && time mpirun -n 2 $0x.x $@ && rm -f $0x.x $0x.cpp; exit
#endif
#ifndef BOOST_MPI3_PROCESS_HPP
#define BOOST_MPI3_PROCESS_HPP

#include "../mpi3/communicator.hpp"

#include <boost/optional.hpp>
using boost::optional;

#include "config/NODISCARD.hpp"

namespace boost{
namespace mpi3{

struct process{
	communicator& comm_;
	int rank_;
	communicator& comm() const{return comm_;}
	int rank() const{return rank_;}
	template<class T>
	optional<T> operator+=(T const& t) &&{
		T val = comm_.reduce_value(t, std::plus<>{}, rank_);
		if(rank_ != comm_.rank()) return {};
		return optional<T>(val);
	}
//	template<class T>
//	std::vector<T> operator|=(T const& t) &&{
//		std::vector<T> ret(comm_.size());
//		comm_.gather_n(&t, 1, ret.begin(), rank_);
	//	comm_.gather_value(t, ret.begin(), rank_);
//		return ret;
//	}
//	template<class T>
//	process&& operator<<(T const& t) &&{
//		comm_.send_value(t, rank_);
//		return std::move(*this);
//	}
	template<class T>
	process&& operator>>(T& t) &&{
		comm_.receive_n(&t, 1, rank_);
	//	comm_.receive_value(t, rank_);
		return std::move(*this);
	}
	template<class T>
	process&& operator&(T& t) &&{
		comm_.broadcast_value(t, rank_);
		return std::move(*this);
	}
};

template<class T>
auto operator<<(process&& p, const T& value) -> decltype(std::move(p << value)){
	return std::move(p << value);
}

template<class T>
auto operator>>(process&& p, T&& value) -> decltype(std::declval<process&>() >> value){
	return p >> value;
}

template<class T> 
process& operator<<(process& self, T const& t){
	self.comm_.send_value(t, self.rank_);
	return self;
}

inline process communicator::operator[](int rank){
	return {*this, rank};
}

template<class T>
auto operator&(communicator& comm, T&& t)
->decltype(comm.all_to_all(begin(std::forward<T>(t))), std::forward<T>(t)){
	assert(t.size() == comm.size());
//	using std::begin;
	auto e = comm.all_to_all(begin(std::forward<T>(t)));
	using std::end;
	assert( e == end(t) );
	return std::forward<T>(t);
}

template<class T> 
NODISCARD("do not ignore result when second argument is const")
auto operator&(communicator& comm, T const& t)
->decltype(comm.all_to_all(t.begin(), std::declval<T>().begin()), T(comm.size())){
	assert(t.size() == comm.size());
	T ret(comm.size()); 
	comm.all_to_all(t.begin(), ret.begin());
	return ret;
}

template<class T>
communicator& operator>>(communicator& comm, T& t){
	comm.receive_n(&t, 1);
//	comm.receive_value(t);
	return comm;
}
template<class T>
std::vector<T> operator|=(communicator& comm, T const& t){
	return comm.all_gather_value(t);
}

template<class T>
std::vector<T> operator|=(process&& self, T const& t){
	return self.comm_.gather_value(t, self.rank_);
}

}}

#ifdef _TEST_BOOST_MPI3_PROCESS

#include "../mpi3/main.hpp"

//#include "alf/boost/multi_array/serialization.hpp"

#include<boost/serialization/vector.hpp>
#include<boost/multi_array.hpp>

#include<boost/archive/xml_iarchive.hpp>
#include<boost/archive/xml_oarchive.hpp>

#include<fstream>


namespace mpi3 = boost::mpi3;
using std::cout;

template<class T>
bool save_xml(std::string file, T const& ma){
	std::ofstream ofs(file.c_str()); assert(ofs);
	boost::archive::xml_oarchive oa(ofs);
	oa << boost::serialization::make_nvp(file.c_str(), ma);
	return ofs.good();
}

template<class T>
bool load_xml(std::string file, T&& ma) try{
	std::ifstream ifs(file.c_str()); assert(ifs);
	boost::archive::xml_iarchive ia(ifs);
	ia >> boost::serialization::make_nvp(file.c_str(), std::forward<T>(ma));
	return ifs.good();
}catch(std::exception& e){
	throw std::runtime_error("cannot load from file `\n" + file + ": line 0:\n', because " + e.what() + ". Are the save and load type compatible?");
}

int mpi3::main(int argc, char* argv[], mpi3::communicator world){
	assert(world.size() == 2);
	if(world.rank() == 0){
		int a = 7;
		world[1] << a;
	}else if(world.rank() == 1){
		int a = -1;
		world >> a; // any source (any tag)
		assert(a == 7);
	}
	switch(world.rank()){
		case 0: {
			bool b = true;
			world[1] << b;
		} break;
		case 1:{
			bool b = false;
			world >> b;
			assert(b == true);
		} break;
	}
	if(world.rank() == 0){
		std::vector<double> v = {
			1, 2, 3,
			4, 5, 6,
			7, 8, 9,
			10, 11, 12
		};
		world[1] << v;
	}else if(world.rank() == 1){
		std::vector<double> v;
		world >> v;
		assert(( v == std::vector<double>{1,2,3,  4, 5, 6,  7,8,9, 10,11,12} ));
	}
	switch(world.rank()){
		case 0:{
			std::vector<bool> v = {false, true, false};
			world[1] << v;
		} break; case 1:{
			std::vector<bool> v;
			world >> v;
			assert(( v == std::vector<bool>{false, true, false} ));
		}
	}
#if 0
	if(world.rank() == 0){
		std::vector<double> v = {1,2,3,  4, 5, 6,  7,8,9,  10,11,12};
		boost::multi_array<double, 2> ma(boost::extents[4][2]);
		std::copy(v.begin(), v.begin() + 8, ma.data());
		save_xml("process.xml", ma);
		world[1] << ma;
	}else if(world.rank() == 1){
		boost::multi_array<double, 2> ma(boost::extents[4][2]);
		world >> ma;
	//	load_xml("process.xml", ma);
		cout << ma.shape()[0] << " " << ma.shape()[1] << '\n';
	//	assert( ma[2][2] == 8. );
	}
#endif
	return 0;
}
#endif
#endif


