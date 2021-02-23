#if COMPILATION_INSTRUCTIONS /* -*- indent-tabs-mode: t -*- */
(echo "#include\""$0"\"" > $0x.cpp) && mpic++ -O3 -std=c++14 -Wfatal-errors -D_TEST_BOOST_MPI3_ERROR_HANDLER $0x.cpp -o $0x.x && time mpirun -n 4 $0x.x $@ && rm -f $0x.cpp; exit
#endif
#ifndef BOOST_MPI3_ERROR_HANDLER_HPP
#define BOOST_MPI3_ERROR_HANDLER_HPP

#define OMPI_SKIP_MPICXX 1  // https://github.com/open-mpi/ompi/issues/5157
#include<mpi.h>

#include "../mpi3/communicator.hpp"

namespace boost{
namespace mpi3{

class fatal{
	operator MPI_Errhandler() const{return MPI_ERRORS_ARE_FATAL;}
};

class code{
	operator MPI_Errhandler() const{return MPI_ERRORS_RETURN;}
};

struct error_handler{
	MPI_Errhandler impl_ = MPI_ERRORS_ARE_FATAL;
	error_handler(MPI_Errhandler impl) : impl_(impl){}
	public:
	error_handler() = default;
//	error_handler(void(*fn)(MPI_Comm*, int* err, ...)){
//		MPI_Comm_create_errhandler(fn, &impl_);
//	}
//	void operator()(communicator& comm, int error) const{comm.call_error_handler(error);}
	~error_handler(){
		if(impl_ != MPI_ERRORS_ARE_FATAL and impl_ != MPI_ERRORS_RETURN)
			MPI_Errhandler_free(&impl_);
	}
	static error_handler const fatal;
	static error_handler const code;
//	static error_handler const exception;
	static void exception(MPI_Comm* /*comm*/, int* err, ...){
		int len = -1;
		char estring[MPI_MAX_ERROR_STRING];
		MPI_Error_string(*err, estring, &len);
		std::string w(estring, estring + len);
		throw std::runtime_error{"error code"};
//		throw boost::mpi3::exception("error code " + std::to_string(*err) + " from comm " + std::to_string(*comm) + ": " + w);
//		throw std::runtime_error("error code " + std::to_string(*err) + " from comm " + std::to_string(*comm) + ": " + w);
	}
};

error_handler const error_handler::fatal(MPI_ERRORS_ARE_FATAL);
error_handler const error_handler::code(MPI_ERRORS_RETURN);

void communicator::set_error_handler(error_handler const& eh){
	int s = MPI_Comm_set_errhandler(impl_, eh.impl_);
	if(s != MPI_SUCCESS) throw std::runtime_error{"cannot set error handler"};
}

error_handler communicator::get_error_handler() const{
	error_handler ret;
	int status = MPI_Comm_get_errhandler(impl_, &ret.impl_);
	if(status != MPI_SUCCESS) throw std::runtime_error("cannot get error handler");
	return ret;
}

}}

#ifdef _TEST_BOOST_MPI3_ERROR_HANDLER

#include "../mpi3/main.hpp"

namespace mpi3 = boost::mpi3;
using std::cout;

int mpi3::main(int, char*[], mpi3::communicator world){

//	error_counter_t ec;
//	mpi3::error_handler<&ehh> ehherr;
//	mpi3::error_handler newerr(eh);
//	world.set_error_handler(ec);
//	world.set_error_handler(newerr);
//	world.set_error_handler(
//	world.error(error::other);
//	cout << ec.count() << '\n';
//	newerr(world, MPI_ERR_OTHER);

//	auto f = world.get_error_handler();

	return 0;
}

#endif
#endif

