#if COMPILATION_INSTRUCTIONS /* -*- indent-tabs-mode: t -*- */
(echo "#include\""$0"\"" > $0x.cpp) && mpic++ -O3 -std=c++14 -Wfatal-errors -Wall -Wextra -D_MAKE_BOOST_SERIALIZATION_HEADER_ONLY `#-lboost_serialization` -D_TEST_BOOST_MPI3_PORT $0x.cpp -o $0x.x && time mpirun -n 4 $0x.x $@ && rm -f $0x.cpp; exit
#endif
#ifndef BOOST_MPI3_PORT_HPP
#define BOOST_MPI3_PORT_HPP

#define OMPI_SKIP_MPICXX 1  // https://github.com/open-mpi/ompi/issues/5157
#include<mpi.h>

#include<stdexcept>
#include<string>

namespace boost{
namespace mpi3{

struct port{
	std::string name_ = ""; // typically this will be tag#0$description#inspiron$port#47425$ifname#172.17.5.240$
	port(){open();}
	port(port const&) = delete;
	port& operator=(port const& other) = delete;
	port(std::string const& name) : name_(name){};// open(name);}
	~port(){ if(is_open()) close(); }
	void open(){
		char name[MPI_MAX_PORT_NAME];
		int status = MPI_Open_port(MPI_INFO_NULL, name);
		name_ = name;
		if(status != 0) throw std::runtime_error("can't open port " + name_);
	}
	void open(std::string const& name){name_ = name;}
	std::string const& name() const{return name_;}
	bool is_open() const{return (name_ != "");}
	void close(){
		int status = MPI_Close_port(name_.c_str());
		if(status != 0) throw std::runtime_error("can't close port" + name_);
		name_ = "";
	}
};

}}

#ifdef _TEST_BOOST_MPI3_PORT

#include "../mpi3/environment.hpp"

using std::cout;
namespace mpi3 = boost::mpi3;

int main(int argc, char* argv[]){
	mpi3::environment env(argc, argv);
	auto world = env.world();
	assert(world.size() > 2);
	switch(world.rank()){
		case 0:{
			mpi3::port p1;
			mpi3::port p2;
			world.send_value(p1.name(), 1);
			world.send_value(p2.name(), 2);
			mpi3::communicator comm1 = env.self().accept(p1, 0);
			mpi3::communicator comm2 = env.self().accept(p2, 0);
			comm1.send_value(1, 0);
			comm2.send_value(2, 0);
		}; break;
		case 1:{
			std::string s;
			world.receive_n(&s, 1, 0);
			mpi3::port p1(s);
			mpi3::communicator comm1 = env.self().connect(p1, 0);
			int data = -1;
			comm1.receive_n(&data, 1, 0);
			assert(data == 1);
		}; break;
		case 2:{
			std::string s;
			world.receive_n(&s, 1);
			mpi3::port p2(s);
			mpi3::communicator comm2 = env.self().connect(p2, 0);
			int data;
			comm2.receive_n(&data, 1, 0); 
			assert(data == 2);
		}; break;
	}
}

#endif
#endif

