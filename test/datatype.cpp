#if COMPILATION_INSTRUCTIONS
mpic++ -g -O3 -Wall -Wextra $0 -o $0x -D_MAKE_BOOST_SERIALIZATION_HEADER_ONLY`#-lboost_serialization`&&mpirun -n 3 valgrind --leak-check=full --show-reachable=yes --error-limit=no                                 --suppressions=communicator_main.cpp.openmpi.supp $0x&&rm $0x;exit
#endif
// © Alfredo A. Correa 2022

#include "../../mpi3/communicator.hpp"
#include "../../mpi3/main.hpp"

#include<complex>
#include<list>
#include<string>

namespace mpi3 = boost::mpi3;

auto mpi3::main(int/*argc*/, char**/*argv*/, mpi3::communicator /*world*/) -> int try{
	using mpi3::detail::is_basic;

	static_assert( is_basic<int>{} );
	static_assert( is_basic<double>{} );
	static_assert( is_basic<mpi3::detail::float_int>{} );

	static_assert( not is_basic<std::string>{} );

	assert( mpi3::detail::basic_datatype<double>{} == MPI_DOUBLE );

	return 0;
}catch(...){
	return 1;
}

