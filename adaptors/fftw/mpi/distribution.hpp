#if COMPILATION_INSTRUCTIONS
#mpicxx -I$HOME/prj/alf $0 -g -o $0x -lfftw3 -lfftw3_mpi &&mpirun -n 4 valgrind $0x;exit
$CXXX $CXXFLAGS -O2 -g `mpicxx -showme:compile|sed 's/-pthread/ /g'` -I$HOME/prj/alf $0 -o $0x `mpicxx -showme:link|sed 's/-pthread/ /g'` -lfftw3 -lfftw3_mpi -lboost_timer&&mpirun -n 4 $0x;exit
#endif

#ifndef MULTI_FFTW_MPI_DISTRIBUTION_HPP
#define MULTI_FFTW_MPI_DISTRIBUTION_HPP

#include <fftw3-mpi.h>

#include<boost/mpi3/communicator.hpp>

#include "../../../array_ref.hpp"

#include <experimental/tuple>

namespace boost{
namespace multi{
namespace fftw{
namespace mpi{

namespace bmpi3 = boost::mpi3;

using difference_type = std::ptrdiff_t;

template<std::ptrdiff_t ElementSize>
class many{
public:
	using difference_type = std::ptrdiff_t;
private:
	difference_type local_count_;
	difference_type local_n0_;
	difference_type local_0_start_;
	static auto sizes(boost::multi::extensions_type_<2> const& ext){
		using std::experimental::apply;
		return apply([](auto... e){return std::array<difference_type, 2>{e.size()...};}, ext);
	}
public:
	many(extensions_type_<2> const& ext, bmpi3::communicator const& comm, difference_type block0 = FFTW_MPI_DEFAULT_BLOCK)
	: local_count_{
		std::max(
			difference_type(
				fftw_mpi_local_size_many(
					2, sizes(ext).data(), ElementSize/sizeof(double),
					block0, comm.get(),
					&local_n0_, &local_0_start_
				)*sizeof(double)/ElementSize
			), 
			difference_type(1)
		)
	}
	{
		static_assert( ElementSize%sizeof(double) == 0 , "!" );
	}
	difference_type local_count() const{return local_count_ + 100;}
	multi::iextension local_extension_0() const{return {local_0_start_, local_0_start_ + local_n0_};}
	multi::iextension local_extension() const{return local_extension_0();}
	bool operator==(many const& other) const{
		return std::tie(this->local_count_, this->local_n0_, this->local_0_start_)
			== std::tie(other.local_count_, other.local_n0_, other.local_0_start_);
	}
	bool operator!=(many const& other) const{return not operator==(other);}
};

template<std::ptrdiff_t ElementSize>
class many_transposed{
public:
	using difference_type = std::ptrdiff_t;
private:
	difference_type local_count_;
	difference_type local_n0_ = 0;
	difference_type local_0_start_ = 0;
	difference_type local_n1_ = 0;
	difference_type local_1_start_ = 0;
	static auto sizes(boost::multi::extensions_type_<2> const& ext){
		using std::experimental::apply;
		return apply([](auto... e){return std::array<difference_type, 2>{e.size()...};}, ext);
	}
public:
	static_assert(ElementSize%sizeof(double)==0, "!");
	many_transposed(
		extensions_type_<2> const& ext, boost::mpi3::communicator const& comm, 
		difference_type block0 = FFTW_MPI_DEFAULT_BLOCK, difference_type block1 = FFTW_MPI_DEFAULT_BLOCK
	) : local_count_{
		std::max(
			difference_type(
				fftw_mpi_local_size_many_transposed(
					2, sizes(ext).data(), ElementSize/sizeof(double),
					block0, block1, comm.get(),
					&local_n0_, &local_0_start_,
					&local_n1_, &local_1_start_
				)*sizeof(double)/ElementSize
			), 
			difference_type(1)
		)
	} {
		static_assert( ElementSize%sizeof(double) == 0 , "!");
		// FFTW_MPI_DEFAULT_BLOCK = (size + comm.size - 1)/comm.size
		assert( local_count() >= local_extension0().size()*local_extension1().size() );
	//	assert( block0*comm.size() >= std::get<0>(ext).size() or block0 == FFTW_MPI_DEFAULT_BLOCK );
	}
	difference_type local_count() const{return local_count_ + 100;}
	multi::iextension local_extension0() const{return {local_0_start_, local_0_start_ + local_n0_};}
	multi::iextension local_extension1() const{return {local_1_start_, local_1_start_ + local_n1_};}
	bool operator==(many_transposed const& other) const{
		return std::tie(this->local_count_, this->local_n0_, this->local_0_start_, this->local_n1_, this->local_1_start_)
			== std::tie(other.local_count_, other.local_n0_, other.local_0_start_, other.local_n1_, other.local_1_start_);
	}
	bool operator!=(many_transposed const& other) const{return not operator==(other);}
};

}}}}

#if not __INCLUDE_LEVEL__

#include<boost/mpi3/main_environment.hpp>
#include<boost/mpi3/ostream.hpp>

#include "../../fftw/mpi/environment.hpp"

namespace bmpi3 = boost::mpi3;
namespace multi = boost::multi;
namespace mpi = multi::fftw::mpi;

int bmpi3::main(int, char*[], mpi3::environment& env){
	multi::fftw::mpi::environment fenv(env);
	auto world = env.world();

	mpi3::ostream os{world};
	
	using std::endl;
	{
		os<< "forced distribution "<<endl;
		mpi::many_transposed<sizeof(double)> dist({12, 43}, world, (12+world.size()-1)/world.size());//533/world.size());
		
		os<< "local element count "<< dist.local_count()            <<endl;
		os<< "local rows          "<< dist.local_extension0().size() <<endl;
		os<< "local extension     "<< dist.local_extension0()        <<endl;
	}
	{
		os<< "automatic distribution "<<std::endl;
		mpi::many_transposed<sizeof(double)> dist({12, 43}, world);//533/world.size());
		
		os<< "local element count "<< dist.local_count()            <<endl;
		os<< "local rows          "<< dist.local_extension0().size() <<endl;
		os<< "local extension     "<< dist.local_extension0()        <<endl;
	}
	mpi::many_transposed<sizeof(double)> forced({12, 43}, world);
	mpi::many_transposed<sizeof(double)> automa({12, 43}, world);
	assert( forced == automa );

	return 0;
}
#endif
#endif

