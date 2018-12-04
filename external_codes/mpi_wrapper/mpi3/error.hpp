#if COMPILATION_INSTRUCTIONS
(echo "#include\""$0"\"" > $0x.cpp) && mpic++ -O3 -std=c++14 -Wall -Wextra -Wfatal-errors -D_TEST_BOOST_MPI3_ERROR $0x.cpp -o $0x.x && time mpirun -n 4 $0x.x $@ && rm -f $0x.cpp; exit
#endif
#ifndef BOOST_MPI3_ERROR_HPP
#define BOOST_MPI3_ERROR_HPP

#define OMPI_SKIP_MPICXX 1  // https://github.com/open-mpi/ompi/issues/5157
#include<mpi.h>

//#include<string>
#include<system_error>

namespace boost{
namespace mpi3{

enum class error : int {//decltype(MPI_SUCCESS) {
	success = MPI_SUCCESS,
	invalid_buffer_pointer = MPI_ERR_BUFFER,
	invalid_count = MPI_ERR_COUNT,
	invalid_datatype = MPI_ERR_TYPE,
	invalid_tag = MPI_ERR_TAG,
	invalid_communicator = MPI_ERR_COMM,
	invalid_rank = MPI_ERR_RANK,
	invalid_root = MPI_ERR_ROOT,
	invalid_group = MPI_ERR_GROUP,
	invalid_operation = MPI_ERR_OP,
	invalid_topology = MPI_ERR_TOPOLOGY,
	illegal_dimension = MPI_ERR_DIMS,
	invalid_dimension = MPI_ERR_DIMS,
	invalid_argument = MPI_ERR_ARG,
	invalid_domain = MPI_ERR_ARG,
	unknown = MPI_ERR_UNKNOWN,
	truncated_message = MPI_ERR_TRUNCATE,
	other = MPI_ERR_OTHER,
	internal = MPI_ERR_INTERN,
	in_status = MPI_ERR_IN_STATUS,
	pending = MPI_ERR_PENDING,
	illegal_request = MPI_ERR_REQUEST,
	last_code = MPI_ERR_LASTCODE
};

struct error_category : std::error_category{
	char const* name() const noexcept override{return "mpi3 wrapper";}
	std::string message(int err) const override{
		char estring[MPI_MAX_ERROR_STRING];
		int len; // int eclass;
	//	MPI_Error_class(err, &eclass);
		MPI_Error_string(err, estring, &len);
		return std::string(estring, len);
	}
	static error_category& instance(){
		static error_category instance;
		return instance;
	}
};

inline std::error_code make_error_code(error err) noexcept{
	return std::error_code(int(err), error_category::instance());
}

#if 0
struct error{
	int impl_;
	error(int code){MPI_Error_class(code, &impl_);}
	int code() const{return impl_;}
	static std::string what(int code){
		int len = -1;
		char estring[MPI_MAX_ERROR_STRING];
		MPI_Error_string(code, estring, &len);
		return std::string(estring, estring + len);
	}
	std::string what() const{return what(impl_);}
	enum code {
		success = MPI_SUCCESS,
		invalid_buffer_pointer = MPI_ERR_BUFFER,
		invalid_count = MPI_ERR_COUNT,
		invalid_datatype = MPI_ERR_TYPE,
		invalid_tag = MPI_ERR_TAG,
		invalid_communicator = MPI_ERR_COMM,
		invalid_rank = MPI_ERR_RANK,
		invalid_root = MPI_ERR_ROOT,
		invalid_group = MPI_ERR_GROUP,
		invalid_operation = MPI_ERR_OP,
		invalid_topology = MPI_ERR_TOPOLOGY,
		illegal_dimension = MPI_ERR_DIMS,
		invalid_dimension = MPI_ERR_DIMS,
		invalid_argument = MPI_ERR_ARG,
		invalid_domain = MPI_ERR_ARG,
		unknown = MPI_ERR_UNKNOWN,
		truncated_message = MPI_ERR_TRUNCATE,
		other = MPI_ERR_OTHER,
		internal = MPI_ERR_INTERN,
		in_status = MPI_ERR_IN_STATUS,
		pending = MPI_ERR_PENDING,
		illegal_request = MPI_ERR_REQUEST,
		last_code = MPI_ERR_LASTCODE
	}; 
};
#endif

}}

namespace std{
	template<> struct is_error_code_enum<::boost::mpi3::error> : true_type{};
}

#ifdef _TEST_BOOST_MPI3_ERROR

#include "../mpi3/main.hpp"
//#include "../mpi3/error_handler.hpp"
#include<ios>

namespace mpi3 = boost::mpi3;
using std::cout;

int mpi3::main(int, char*[], mpi3::communicator world){

	std::error_code ec = mpi3::error::invalid_buffer_pointer;
	assert( ec == mpi3::error::invalid_buffer_pointer);
	assert( ec != std::io_errc::stream );	
	
	try{
		world.broadcast_n((int*)nullptr, 0, -1);
	}catch(std::system_error const& e){
		cout << e.code() <<'\n';
		cout << e.code().message() <<'\n';
		cout << e.what() <<'\n';
	}

	return 0;

}

#endif
#endif

