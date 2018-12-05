#ifndef BOOST_MPI3_DETAIL_CALL_HPP
#define BOOST_MPI3_DETAIL_CALL_HPP

#define OMPI_SKIP_MPICXX 1  // https://github.com/open-mpi/ompi/issues/5157
#include<mpi.h> // MPI_MAX_PROCESSOR_NAME

#include<string>

namespace boost{
namespace mpi3{
namespace detail{

template<int(*F)(char*, int*)> std::string call(){
	int len = -1;
	char name[MPI_MAX_PROCESSOR_NAME];
	F(name, &len);
	return std::string(name, name + len);
}

}}}

#endif

