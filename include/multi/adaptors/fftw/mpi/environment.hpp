#if COMPILATION_INSTRUCTIONS
#mpicxx -I$HOME/prj/alf $0 -g -o $0x -lfftw3 -lfftw3_mpi &&mpirun -n 2 valgrind $0x;exit
$CXXX $CXXFLAGS -O2 -g `mpicxx -showme:compile|sed 's/-pthread/ /g'` -I$HOME/prj/alf $0 -o $0x `mpicxx -showme:link|sed 's/-pthread/ /g'` -lfftw3 -lfftw3_mpi -lboost_timer&&mpirun -n 2 $0x;exit
#endif

#ifndef MULTI_FFTW_MPI_ENVIRONMENT_HPP
#define MULTI_FFTW_MPI_ENVIRONMENT_HPP

#include <fftw3-mpi.h>

#include<boost/mpi3/communicator.hpp>

#include "../../../array_ref.hpp"

#include <experimental/tuple>

namespace boost{
namespace multi{
namespace fftw{
namespace mpi{

namespace bmpi3 = boost::mpi3;

struct environment{
	explicit environment(bmpi3::environment&){fftw_mpi_init();}
	~environment(){fftw_mpi_cleanup();}
};

}}}}

#if not __INCLUDE_LEVEL__

#include<boost/mpi3/main_environment.hpp>

namespace bmpi3 = boost::mpi3;
namespace multi = boost::multi;

int bmpi3::main(int, char*[], mpi3::environment& env){
	multi::fftw::mpi::environment fenv(env);
	return 0;
}
#endif
#endif

