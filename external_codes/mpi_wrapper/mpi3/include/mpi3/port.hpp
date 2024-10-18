// Copyright 2018-2023 Alfredo A. Correa

#ifndef BMPI3_PORT_HPP
#define BMPI3_PORT_HPP

#include<mpi.h>

#include<stdexcept>
#include<string>

namespace boost {
namespace mpi3 {

#if not defined(EXAMPI)
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
		name_.resize(MPI_MAX_PORT_NAME, '\0');
		MPI_(Open_port)(MPI_INFO_NULL, name_.data());
	}
	void open(std::string const& name) {name_ = name;}

	std::string const& name() const{return name_;}

	bool is_open() const {return not name_.empty();}
	void close() {
		MPI_(Close_port)(name_.c_str());
		name_ = "";
	}
};
#endif

}  // end namespace mpi3
}  // end namespace boost
#endif
