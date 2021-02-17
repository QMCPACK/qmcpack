#if COMPILATION// -*- indent-tabs-mode: t -*-
OMPI_CXX=$CXX mpic++ $0 -o $0x&&mpirun -n 4 $0x&&rm $0x;exit
#endif
// © Alfredo A. Correa 2017-2020

#ifndef BOOST_MPI3_ERROR_HPP
#define BOOST_MPI3_ERROR_HPP

#define OMPI_SKIP_MPICXX 1  // https://github.com/open-mpi/ompi/issues/5157
#include<mpi.h>

#include<system_error>

namespace boost{
namespace mpi3{

static_assert(sizeof(MPI_SUCCESS) <= sizeof(int), "!");

enum class error : int {//decltype(MPI_SUCCESS){
	success                = MPI_SUCCESS,
	invalid_buffer_pointer = MPI_ERR_BUFFER,
	invalid_count          = MPI_ERR_COUNT,
	invalid_datatype       = MPI_ERR_TYPE,
	invalid_tag            = MPI_ERR_TAG,
	invalid_communicator   = MPI_ERR_COMM,
	invalid_rank           = MPI_ERR_RANK,
	invalid_root           = MPI_ERR_ROOT,
	invalid_group          = MPI_ERR_GROUP,
	invalid_operation      = MPI_ERR_OP,
	invalid_topology       = MPI_ERR_TOPOLOGY,
	illegal_dimension      = MPI_ERR_DIMS,
	invalid_dimension      = MPI_ERR_DIMS,
	invalid_argument       = MPI_ERR_ARG,
	invalid_domain         = MPI_ERR_ARG,
	unknown                = MPI_ERR_UNKNOWN,
	truncated_message      = MPI_ERR_TRUNCATE,
	other                  = MPI_ERR_OTHER,
	internal               = MPI_ERR_INTERN,
	in_status              = MPI_ERR_IN_STATUS,
	pending                = MPI_ERR_PENDING,
	illegal_request        = MPI_ERR_REQUEST,
	last_code              = MPI_ERR_LASTCODE
};

auto inline string(enum error err){
	char estring[MPI_MAX_ERROR_STRING];
	int len;
	MPI_Error_string(static_cast<int>(err), estring, &len);
	return std::string(estring, len);
}

struct error_category : std::error_category{
	char const* name() const noexcept override{return "mpi3 wrapper";}
	std::string message(int err) const override{return string(static_cast<enum error>(err));}
	static error_category& instance(){
		static error_category instance;
		return instance;
	}
};

inline std::error_code make_error_code(error err) noexcept{
	return std::error_code(int(err), error_category::instance());
}

}}

namespace std{
	template<> struct is_error_code_enum<::boost::mpi3::error> : true_type{};
}

#if not __INCLUDE_LEVEL__ // def _TEST_BOOST_MPI3_ERROR

#include "../mpi3/main.hpp"

namespace mpi3 = boost::mpi3;
using std::cout;

int mpi3::main(int, char*[], mpi3::communicator world){

	std::error_code ec = mpi3::error::invalid_buffer_pointer;
	assert( ec == mpi3::error::invalid_buffer_pointer);
	assert( ec != std::io_errc::stream );	
	
	try{
		world.broadcast_n((int*)nullptr, 0, -1);
	}catch(std::system_error const& e){
		cout
			<<"code: "   << e.code()           <<'\n'
			<<"message: "<< e.code().message() <<'\n'
			<<"what: "   << e.what()           <<'\n'
		;
	}

	return 0;

}

#endif
#endif

