#ifndef BOOST_MPI3_EXCEPTION
#define BOOST_MPI3_EXCEPTION

#include<stdexcept>

namespace boost{
namespace mpi3{

struct exception : std::runtime_error{
	using runtime_error::runtime_error;
//	std::string what_;
//public:
/*	exception(const char* what) : what_(what){}
	virtual const char* what() const noexcept{
		return what_.c_str();
	}*/
	virtual ~exception(){}
};

struct invalid_communicator : exception{
	using exception::exception;
	virtual ~invalid_communicator(){}
};

}}

#endif

