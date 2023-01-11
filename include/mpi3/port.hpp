/* -*- indent-tabs-mode: t -*- */

#ifndef BOOST_MPI3_PORT_HPP
#define BOOST_MPI3_PORT_HPP

// #define OMPI_SKIP_MPICXX 1  // https://github.com/open-mpi/ompi/issues/5157
#include<mpi.h>

#include<stdexcept>
#include<string>

namespace boost {
namespace mpi3 {

struct port {
	// NOLINTNEXTLINE(misc-non-private-member-variables-in-classes) TODO(correaa)
	std::string name_;  // typically this will be something like tag#0$description#inspiron$port#47425$ifname#172.17.5.240$

	port() {open();}
	port(port const&) = delete;
	port(port     &&) = delete;

	port& operator=(port const&) = delete;
	port& operator=(port     &&) = delete;

	explicit port(std::string name) : name_{std::move(name)} {};

	~port() noexcept{ try{if(is_open()) {close();}}catch(...){} }

	void open() {
		std::array<char, MPI_MAX_PORT_NAME> name_buffer{};
		MPI_(Open_port)(MPI_INFO_NULL, name_buffer.data());
		name_ = std::string{name_buffer.data()};
	}
	void open(std::string const& name) {name_ = name;}

	std::string const& name() const{return name_;}

	bool is_open() const {return not name_.empty();}
	void close() {
		MPI_(Close_port)(name_.c_str());
		name_ = "";
	}
};

}  // end namespace mpi3
}  // end namespace boost

//#ifdef _TEST_BOOST_MPI3_PORT

//#include "../mpi3/environment.hpp"

//using std::cout;
//namespace mpi3 = boost::mpi3;

//int main(int argc, char* argv[]){
//	mpi3::environment env(argc, argv);
//	auto world = env.world();
//	assert(world.size() > 2);
//	switch(world.rank()){
//		break; case 0:{
//			mpi3::port p1;
//			mpi3::port p2;
//			world.send_value(p1.name(), 1);
//			world.send_value(p2.name(), 2);
//			mpi3::communicator comm1 = env.self().accept(p1, 0);
//			mpi3::communicator comm2 = env.self().accept(p2, 0);
//			comm1.send_value(1, 0);
//			comm2.send_value(2, 0);
//		};
//		break; case 1:{
//			std::string s;
//			world.receive_n(&s, 1, 0);
//			mpi3::port p1(s);
//			mpi3::communicator comm1 = env.self().connect(p1, 0);
//			int data = -1;
//			comm1.receive_n(&data, 1, 0);
//			assert(data == 1);
//		};
//		break; case 2:{
//			std::string s;
//			world.receive_n(&s, 1);
//			mpi3::port p2(s);
//			mpi3::communicator comm2 = env.self().connect(p2, 0);
//			int data;
//			comm2.receive_n(&data, 1, 0); 
//			assert(data == 2);
//		};
//	}
//}

//#endif
#endif

