#if COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;-*-
 ln -sf $0 $0.cpp;mpicxx -g -I$HOME/prj/alf $0.cpp -o $0x -lfftw3 -lfftw3_mpi&&time mpirun -n 4 $0x&&rm $0x $0.cpp;exit
#ln -sf $0 $0.cpp;mpicxx -g -I$HOME/prj/alf $0.cpp -o $0x -lfftw3 -lfftw3_mpi&&time mpirun -n 4 valgrind --leak-check=full --track-origins=yes --show-leak-kinds=all --suppressions=$HOME/prj/alf/boost/mpi3/test/communicator_main.cpp.openmpi.supp --error-exitcode=1 $0x&&rm $0x $0.cpp;exit
#endif
// © Alfredo A. Correa 2020
// apt-get install libfftw3-mpi-dev
// compile with: mpicc simple_mpi_example.c  -Wl,-rpath=/usr/local/lib -lfftw3_mpi -lfftw3 -o simple_mpi_example */

#include "../../array.hpp"
#include "../../config/NODISCARD.hpp"

#include<boost/mpi3/communicator.hpp>
#include<boost/mpi3/environment.hpp>

#include "../fftw.hpp"

#include <fftw3-mpi.h>

namespace boost{
namespace multi{
namespace fftw{

template<typename T>
struct allocator : std::allocator<T>{
	template <typename U> struct rebind{using other = fftw::allocator<U>;};
	NODISCARD("to avoid memory leak") 
	T* allocate(std::size_t n){ return static_cast<T*>(fftw_malloc(sizeof(T)*n));}
	void deallocate(T* data, std::size_t){fftw_free(data);}
};

namespace mpi{

struct environment{
	 environment(){fftw_mpi_init();}
	~environment(){fftw_mpi_cleanup();}
};

template<class T, multi::dimensionality_type D, class Alloc = fftw::allocator<T>> 
struct array;

namespace bmpi3 = boost::mpi3;

template<class T, class Alloc>
struct array<T, multi::dimensionality_type{2}, Alloc>{
	using element_type = T;

	mutable bmpi3::communicator comm_;
	Alloc           alloc_;

	typename std::allocator_traits<Alloc>::size_type                local_count_;
	array_ptr<T, 2, typename std::allocator_traits<Alloc>::pointer> local_ptr_;
	ptrdiff_t                                                       n0_;

	static std::pair<typename std::allocator_traits<Alloc>::size_type, multi::extensions_type_<2>> 
	local_2d(multi::extensions_type_<2> ext, boost::mpi3::communicator const& comm){
		ptrdiff_t local_n0, local_0_start;
		auto count = fftw_mpi_local_size_2d(std::get<0>(ext).size(), std::get<1>(ext).size(), comm.get(), &local_n0, &local_0_start);
		assert( count >= local_n0*std::get<1>(ext).size() );
		return {count, {{local_0_start, local_0_start + local_n0}, std::get<1>(ext)}};
	}
	static auto local_count_2d(multi::extensions_type_<2> ext, boost::mpi3::communicator const& comm){
		return local_2d(ext, comm).first;
	}
	static auto local_extension_2d(multi::extensions_type_<2> ext, boost::mpi3::communicator const& comm){
		return local_2d(ext, comm).second;
	}
	array(multi::extensions_type_<2> ext, bmpi3::communicator comm = mpi3::environment::self(), Alloc alloc = {}) :
		comm_{std::move(comm)},
		alloc_{alloc},
		local_count_{local_count_2d(ext, comm_)},
		local_ptr_  {alloc_.allocate(local_count_), local_extension_2d(ext, comm_)},
		n0_{multi::layout_t<2>(ext).size()}
	{
		if(not std::is_trivially_default_constructible<element_type>{})
			adl_alloc_uninitialized_default_construct_n(alloc_, local_ptr_->base(), local_ptr_->num_elements());
	}
	bmpi3::communicator& comm() const&{return comm_;}
	array(array const& other) :
		comm_       {other.comm_},
		alloc_      {other.alloc_},
		local_count_{other.local_count_},
		local_ptr_  {alloc_.allocate(local_count_), local_extension_2d(other.extensions(), comm_)},
		n0_{multi::layout_t<2>(other.extensions()).size()}
	{
		local_cutout() = other.local_cutout();
	}
	array(array&& other) :
		comm_       {std::move(other.comm_)},
		alloc_      {std::move(other.alloc_)},
		local_count_{std::exchange(other.local_count_, 0)},
		local_ptr_  {std::exchange(other.local_ptr_, nullptr)},
		n0_{multi::layout_t<2>(other.extensions()).size()}
	{}
	explicit array(multi::array<T, 2> const& other, bmpi3::communicator comm = mpi3::environment::self(), Alloc alloc = {}) :
		array(other.extensions(), comm, alloc)
	{
		local_cutout() = other.stenciled(std::get<0>(local_cutout().extensions()), std::get<1>(local_cutout().extensions()));
	}
	bool empty() const{return extensions().num_elements();}
	array_ref <T, 2> local_cutout()      &{return *local_ptr_;}
	array_cref<T, 2> local_cutout() const&{return *local_ptr_;}
	ptrdiff_t        local_count() const&{return  local_count_;}
	multi::extensions_type_<2> extensions() const&{return {n0_, std::get<1>(local_cutout().extensions())};}
	ptrdiff_t num_elements() const&{return multi::layout_t<2>(extensions()).num_elements();}
	operator multi::array<T, 2>() const&{ static_assert( std::is_trivially_copy_assignable<T>{}, "!" );
		multi::array<T, 2> ret(extensions(), alloc_);
		comm_.all_gatherv_n(local_cutout().data_elements(), local_cutout().num_elements(), ret.data_elements());
		return ret;
	}
	array& operator=(multi::array<T, 2> const& other) &{
		if(other.extensions() == extensions()) local_cutout() = other.stenciled(std::get<0>(local_cutout().extensions()), std::get<1>(local_cutout().extensions()));
		else{
			array tmp{other};
			std::swap(*this, tmp);
		}
		return *this;
	}
	bool operator==(multi::array<T, 2> const& other) const&{
		if(other.extensions() != extensions()) return false;
		return comm_&=(local_cutout() == other.stenciled(std::get<0>(local_cutout().extensions()), std::get<1>(local_cutout().extensions())));
	}
	friend bool operator==(multi::array<T, 2> const& other, array const& self){
		return self.operator==(other);
	}
	bool operator==(array<T, 2> const& other) const&{assert(comm_==other.comm_);
		return comm_&=(local_cutout() == other.local_cutout());
	}
	array& operator=(array const& other)&{
		if(other.extensions() == this->extensions() and other.comm_ == other.comm_)
			local_cutout() = other.local_cutout();
		else assert(0);
		return *this;
	}
	~array() noexcept{alloc_.deallocate(local_cutout().data_elements(), local_count_);}
};

array<std::complex<double>, 2>& dft(array<std::complex<double>, 2> const& A, array<std::complex<double>, 2>& B, fftw::sign s){
	assert( A.extensions() == B.extensions() );
	assert( A.comm() == B.comm() );
	fftw_plan p = fftw_mpi_plan_dft_2d(
		std::get<0>(A.extensions()).size(), std::get<1>(A.extensions()).size(), 
		(fftw_complex *)A.local_cutout().data_elements(), (fftw_complex *)B.local_cutout().data_elements(), 
		A.comm().get(),
		s, FFTW_ESTIMATE
	);
	fftw_execute(p);
	fftw_destroy_plan(p);
	return B;
}

array<std::complex<double>, 2>& dft_forward(array<std::complex<double>, 2> const& A, array<std::complex<double>, 2>& B){
	return dft(A, B, fftw::forward);
}

array<std::complex<double>, 2> dft_forward(array<std::complex<double>,2> const& A){
	array<std::complex<double>, 2> ret(A.extensions()); dft_forward(A, ret); return ret;
}

}}}}

#if not __INCLUDE_LEVEL__

#include<boost/mpi3/main.hpp>
#include<boost/mpi3/environment.hpp>
#include<boost/mpi3/ostream.hpp>
#include "../fftw.hpp"

namespace mpi3 = boost::mpi3;
namespace multi = boost::multi;

int mpi3::main(int, char*[], mpi3::communicator world){
	multi::fftw::mpi::environment fenv;

	multi::fftw::mpi::array<std::complex<double>, 2> A({41, 321}, world);

	mpi3::ostream os{world};
	os<< "global sizes" << std::get<0>(A.extensions()) <<'x'<< std::get<1>(A.extensions()) <<' '<< A.num_elements() <<std::endl;
	os<< A.local_cutout().extension() <<'x'<< std::get<1>(A.local_cutout().extensions()) <<"\t#="<< A.local_cutout().num_elements() <<" allocated "<< A.local_count() <<std::endl;

	{
		auto x = A.local_cutout().extensions();
		for(auto i : std::get<0>(x))
			for(auto j : std::get<1>(x))
				A.local_cutout()[i][j] = std::complex<double>(i + j, i + 2*j + 1)/std::abs(std::complex<double>(i + j, i + 2*j + 1));
	}
	
	multi::array<std::complex<double>, 2> A2 = A;
	assert( A2 == A );
	
	using multi::fftw::dft_forward;

	dft_forward(A , A );
	dft_forward(A2, A2);

	{
		auto x = A.local_cutout().extensions();
		for(auto i : std::get<0>(x))
			for(auto j : std::get<1>(x))
				if(not( std::abs(A.local_cutout()[i][j] - A2[i][j]) < 1e-12 )){
					std::cout << A.local_cutout()[i][j] - A2[i][j] <<' '<< std::abs(A.local_cutout()[i][j] - A2[i][j]) << std::endl;
				}
	}
	return 0;
}
#endif


