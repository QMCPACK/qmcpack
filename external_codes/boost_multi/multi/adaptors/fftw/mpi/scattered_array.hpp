#if COMPILATION// -*- indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil; -*-
$CXXX $CXXFLAGS `mpicxx -showme:compile|sed 's/-pthread/ /g'` -I$HOME/prj/alf $0 -o $0x `mpicxx -showme:link|sed 's/-pthread/ /g'` -lfftw3 -lfftw3_mpi&&mpirun -n 4 $0x&&rm $0x;exit
#ln -sf $0 $0.cpp;mpicxx -g -I$HOME/prj/alf $0.cpp -o $0x -lfftw3 -lfftw3_mpi&&mpirun -n 4 valgrind --leak-check=full --track-origins=yes --show-leak-kinds=all --suppressions=$HOME/prj/alf/boost/mpi3/test/communicator_main.cpp.openmpi.supp --error-exitcode=1 $0x&&rm $0x $0.cpp;exit
#endif
// © Alfredo A. Correa 2020
// apt-get install libfftw3-mpi-dev
// compile with: mpicc simple_mpi_example.c  -Wl,-rpath=/usr/local/lib -lfftw3_mpi -lfftw3 -o simple_mpi_example */

#ifndef MULTI_ADAPTOR_FFTW_MPI_SCATTERED_ARRAY_HPP
#define MULTI_ADAPTOR_FFTW_MPI_SCATTERED_ARRAY_HPP

#include "../mpi/distribution.hpp"
#include "boost/mpi3/process.hpp"

namespace boost{
namespace multi{
namespace fftw{
namespace mpi{

namespace bmpi3 = boost::mpi3;

template<class T, multi::dimensionality_type D, class Alloc = std::allocator<T>> // cannot use fftw::allocator<T> as default because it produces error in nvcc:   `template<class _Tp> using __pointer = typename _Tp::pointer’ is protected within this context`
class scattered_array;

template<class T, multi::dimensionality_type D, class Alloc = std::allocator<T>> // cannot use fftw::allocator<T> as default because it produces error in nvcc:   `template<class _Tp> using __pointer = typename _Tp::pointer’ is protected within this context`
class gathered_array;

template<class T, class Alloc>
struct array{
	using local_distrubution_type = many_transposed<sizeof(T)>;
	using local_allocator_type = Alloc;
	using local_pointer_type = typename std::allocator_traits<local_allocator_type>::pointer;
protected:
	local_distrubution_type local_distrubution_;
	local_allocator_type alloc_;
	local_pointer_type local_data_;
	multi::iextension first_ext_;
	multi::iextension second_ext_;
public:
	array(
		multi::extensions_type_<2> exts, bmpi3::communicator comm, 
		difference_type block0 = FFTW_MPI_DEFAULT_BLOCK, difference_type block1 = FFTW_MPI_DEFAULT_BLOCK,
		Alloc alloc = {}
	) :	
		local_distrubution_{exts, comm, block0, block1},
		alloc_{alloc}, 
		local_data_{alloc_.allocate(local_distrubution_.local_count())},
		first_ext_{std::get<0>(exts)},
		second_ext_{std::get<1>(exts)}
	{}
	~array() noexcept{alloc_.deallocate(local_data_, local_distrubution_.local_count());}
	auto local_cutout()      &{return array_ref <T, 2, local_pointer_type>(local_data_, local_distrubution_.local_extension0()*local_distrubution_.local_extension1());}//.rotated();}
	auto local_cutout() const&{return array_cref<T, 2, local_pointer_type>(local_data_, local_distrubution_.local_extension0()*local_distrubution_.local_extension1());}//.rotated();}
};

template<class T, class Alloc>
class gathered_array<T, 2, Alloc> : public array<T, Alloc>{
	bmpi3::communicator comm_;
public:
	gathered_array(multi::extensions_type_<2> exts, bmpi3::communicator comm, Alloc alloc = {}) :
		array<T, Alloc>{exts, comm, std::get<0>(exts).size(), std::get<1>(exts).size(), alloc}, 
		comm_{std::move(comm)}
	{}
	scattered_array<T, 2, Alloc> scatter() const{
		scattered_array<T, 2, Alloc> other({this->first_ext_, this->second_ext_}, comm_);
		auto p = fftw_mpi_plan_many_transpose(
			this->second_ext_.size(), this->first_ext_.size(), 
			sizeof(T)/sizeof(double), 
			this->second_ext_.size(), FFTW_MPI_DEFAULT_BLOCK, 
			reinterpret_cast<double*>(const_cast<T*>(this->local_cutout().base())), 
			reinterpret_cast<double*>(               other.local_cutout().base() ),
			comm_.get(), FFTW_ESTIMATE | FFTW_MPI_TRANSPOSED_IN
		);
		fftw_execute(p);
		fftw_destroy_plan(p);
		return other;
	}
};

template<class T, class Alloc>
class scattered_array<T, 2, Alloc>{
public:
	using local_distrubution_type = many<sizeof(T)>;
	using local_allocator_type = Alloc;
	using local_pointer_type = typename std::allocator_traits<local_allocator_type>::pointer;
private:
	local_distrubution_type local_distribution_;
	local_allocator_type alloc_;
	local_pointer_type local_data_;
	multi::iextension first_ext_;
	multi::iextension second_ext_;
	mutable bmpi3::communicator comm_;
public:
	scattered_array(multi::extensions_type_<2> exts, bmpi3::communicator comm, Alloc alloc = {}) : 
		local_distribution_{exts, comm},
		alloc_{alloc},
		local_data_{alloc_.allocate(local_distribution_.local_count())},
		first_ext_{std::get<0>(exts)},
		second_ext_{std::get<1>(exts)},
		comm_{std::move(comm)}
	{}
	~scattered_array() noexcept{alloc_.deallocate(local_data_, local_distribution_.local_count());}

	array_ref <T, 2, local_pointer_type> local_cutout()      &{return array_ref <T, 2, local_pointer_type>(local_data_, local_distribution_.local_extension_0()*second_ext_);}
	array_cref<T, 2, local_pointer_type> local_cutout() const&{return array_cref<T, 2, local_pointer_type>(local_data_, local_distribution_.local_extension_0()*second_ext_);}

	mpi::gathered_array<T, 2> gather() const{
		mpi::gathered_array<T, 2> other({first_ext_, second_ext_}, comm_);
		auto p = fftw_mpi_plan_many_transpose(
			first_ext_.size(), second_ext_.size(),
		//	std::get<0>(this->extensions()).size(), std::get<1>(this->extensions()).size(), 
			sizeof(T)/sizeof(double), 
			FFTW_MPI_DEFAULT_BLOCK, second_ext_.size(), //this->size(),
			reinterpret_cast<double*>(const_cast<T*>(this->local_cutout().base())), 
			reinterpret_cast<double*>(               other.local_cutout().base() ),
			comm_.get(), FFTW_ESTIMATE | FFTW_MPI_TRANSPOSED_OUT
		);
		fftw_execute(p);
		fftw_destroy_plan(p);
		return other;
	}
	bool operator==(scattered_array const& other) const{return comm_&=(local_cutout() == other.local_cutout());}
	bool operator!=(scattered_array const& other) const{return not operator==(other);}	
};

}}}}

#if 0
template<class T, class Alloc>
class scattered_array<T, multi::dimensionality_type{2}, Alloc>{
public:
	using local_allocator_type = Alloc;
	using local_pointer_t      = typename std::allocator_traits<local_allocator_type>::pointer;
private:
	using local_distrubution_type = distribution<sizeof(T)>;
	using layout_type = layout_t<T, 2>;
	Alloc           alloc_ ;
	local_pointer_t local_data_; //	typename boost::multi::array_ptr<T, 2, local_pointer_t>   local_ptr_;
public:
	scattered_array(multi::extensions_type_<2> ext, bmpi3::communicator comm = mpi3::environment::self(), Alloc alloc = {}) :
		layout_t<T, 2>(ext, comm),
		alloc_{alloc},
		local_data_{std::allocator_traits<Alloc>::allocate(alloc_, scattered_array::local_count())}//,
	{
		if(not std::is_trivially_default_constructible<typename scattered_array::element_type>{})
			adl_alloc_uninitialized_default_construct_n(alloc_, local_cutout().data_elements()/*local_ptr_->base()*/, local_cutout().num_elements());//local_ptr_->num_elements());
	}
	scattered_array(scattered_array const& other) :
		layout_t<T, 2> {other},
		alloc_      {other.alloc_},
		local_data_ {std::allocator_traits<Alloc>::allocate(alloc_, layout_type::local_count())}
	{
		scoped_barrier(other.comm());
		local_cutout() = other.local_cutout();
	/*
		auto p1 = fftw_mpi_plan_many_transpose(
			std::get<0>(this->extensions()).size(), std::get<1>(this->extensions()).size(), sizeof(T)/sizeof(double), 
			other.block(), this->block(),
			reinterpret_cast<double*>(const_cast<T*>(other.local_cutout().data_elements())), 
			reinterpret_cast<double*>(               this->local_cutout().data_elements() ),
			this->comm().get(), FFTW_ESTIMATE
		);
		auto p2 = fftw_mpi_plan_many_transpose(
			std::get<1>(this->extensions()).size(), std::get<0>(this->extensions()).size(), sizeof(T)/sizeof(double), 
			other.block(), this->block(),
			reinterpret_cast<double*>(               this->local_cutout().data_elements()), 
			reinterpret_cast<double*>(               this->local_cutout().data_elements()),
			this->comm().get(), FFTW_ESTIMATE
		);
		fftw_execute(p1);
		fftw_execute(p2);
		fftw_destroy_plan(p2);
		fftw_destroy_plan(p1);
	*/
	}
	scattered_array(scattered_array&& other) : // intel calls this function to return from a function
		layout_type{std::exchange(static_cast<layout_type&>(other), layout_type(multi::extensions_type_<2>{}, other.comm()))},
		alloc_     {std::move(other.alloc_)},
		local_data_{other.local_data_}
	{
		assert(not other.extensions());
		assert(other.local_count() == 0 );
	}
	
	friend std::ostream& operator<<(std::ostream& os, scattered_array const& self){
		for(int r = 0; r != self.comm().size(); ++r){
			if(self.comm().rank() == r){
				if(auto x = self.local_cutout().extensions())
					for(auto i : std::get<0>(x)){
						for(auto j : std::get<1>(x))
							os<< self.local_cutout()[i][j] <<' ';
						os<<std::endl;
					}
			}
			self.comm().barrier();
		}
		return os;
	}

	array_ref <T, 2, local_pointer_t> local_cutout()      &//{return *local_ptr_;}
		{return array_ref <T, 2, local_pointer_t>(local_data_, this->local_extensions());}
	array_cref<T, 2, local_pointer_t> local_cutout() const&//{return *local_ptr_;}
		{return array_cref<T, 2, local_pointer_t>(local_data_, this->local_extensions());}

	local_pointer_t local_data(){return local_data_;}
	typename std::pointer_traits<local_pointer_t>::template rebind<T const> local_data() const{return local_data_;}

	auto extensions() const{return this->global_extensions();}

	operator multi::array<T, 2>() const&{ 
		static_assert( std::is_trivially_copy_assignable<T>{}, "!" );
		multi::array<T, 2> ret(this->global_extensions(), 1., alloc_);
		this->comm().all_gatherv_n(local_data_, local_cutout().num_elements(), ret.data_elements());
		return ret;
	}
	
	mpi::gathered_array<T, 2> gather() const{
		mpi::gathered_array<T, 2> other(this->extensions(), this->comm());
		this->comm_.gatherv_n(local_cutout().data_elements(), local_cutout().num_elements(), other.data_elements());
		static_assert( std::is_trivially_copy_assignable<T>{} and sizeof(T)%sizeof(double)==0, "!");

	/*	{
			fftw_plan p = fftw_mpi_plan_many_transpose(
				std::get<0>(this->extensions()).size(), std::get<1>(this->extensions()).size(), sizeof(T)/sizeof(double), 
				this->block(), std::get<0>(this->extensions()).size(),
				reinterpret_cast<double*>(const_cast<T*>(local_cutout().data_elements())), 
				reinterpret_cast<double*>(ret.data_elements()),
				this->comm().get(), FFTW_ESTIMATE | FFTW_MPI_TRANSPOSED_IN | FFTW_MPI_TRANSPOSED_OUT
			);
			fftw_execute(p);
			fftw_destroy_plan(p);
		}*/
	
		auto p1 = fftw_mpi_plan_many_transpose(
			std::get<0>(this->extensions()).size(), std::get<1>(this->extensions()).size(), 
			sizeof(T)/sizeof(double), 
			FFTW_MPI_DEFAULT_BLOCK, this->size(),
			reinterpret_cast<double*>(const_cast<T*>(this->local_cutout().data_elements())), 
			reinterpret_cast<double*>(               other.data_elements() ),
			this->comm().get(), FFTW_ESTIMATE
		);

		auto p2 = fftw_mpi_plan_many_transpose(
			std::get<1>(this->extensions()).size(), std::get<0>(this->extensions()).size(), 
			sizeof(T)/sizeof(double), 
			other.block(), other.block(),
			reinterpret_cast<double*>(               other.data_elements()), 
			reinterpret_cast<double*>(               other.data_elements()),
			this->comm().get(), FFTW_ESTIMATE
		);
		fftw_execute(p1);
		fftw_execute(p2);
		fftw_destroy_plan(p2);
		fftw_destroy_plan(p1);

		return other;
	}

	explicit scattered_array(multi::array<T, 2> const& other, bmpi3::communicator comm = mpi3::environment::self(), Alloc alloc = {}) :
		scattered_array(other.extensions(), comm, alloc)
	{
		local_cutout() = other.stenciled(std::get<0>(local_cutout().extensions()), std::get<1>(local_cutout().extensions()));
	}
//	bool operator==(array<T, 2> const& other) const&{assert(comm()==other.comm());
//		return comm()&=(local_cutout() == other.local_cutout());
//	}
//	bool operator!=(array<T, 2> const& other) const&{return not(*this==other);}
	ptrdiff_t num_elements() const&{return multi::layout_t<2>(extensions()).num_elements();}
	layout_type layout() const{return *this;}
	~scattered_array() noexcept{if(this->local_count()) alloc_.deallocate(local_data_, this->local_count());}
	
	scattered_array& operator=(scattered_array const& other)&{
		assert(this->comm() == other.comm());
		if(this->extensions() == other.extensions()){
			fftw_plan p = fftw_mpi_plan_many_transpose(
				std::get<0>(this->extensions()).size(), std::get<1>(this->extensions()).size(), sizeof(T)/sizeof(double), 
				other.block(), this->block(),
				reinterpret_cast<double*>(const_cast<T*>(other.local_cutout().data_elements())), 
				reinterpret_cast<double*>(               this->local_cutout().data_elements() ),
				this->comm().get(), FFTW_ESTIMATE
			);
			fftw_execute(p);
			fftw_destroy_plan(p);
		}else assert(0);
		return *this;
	}
#if 0
private:
	typename std::allocator_traits<Alloc>::size_type 
	     local_count_2d    (multi::extensions_type_<2> ext){return local_2d(ext).first; }
	auto local_extension_2d(multi::extensions_type_<2> ext){return local_2d(ext).second;}
public:
	Alloc get_allocator() const{return alloc_;}
	array(bmpi3::communicator comm = mpi3::environment::self(), Alloc alloc = {}) :
		comm_{std::move(comm)},
		alloc_{alloc},
		local_count_{local_count_2d(multi::extensions_type_<2>{})},
		local_ptr_  {alloc_.allocate(local_count_), local_extension_2d(multi::extensions_type_<2>{})},
		n0_{multi::layout_t<2>(multi::extensions_type_<2>{}).size()}
	{}
	bool empty() const{return extensions().num_elements();}
	array_ref <T, 2> local_cutout()      &{return *local_ptr_;}
	array_cref<T, 2> local_cutout() const&{return *local_ptr_;}
	ptrdiff_t        local_count() const&{return  local_count_;}
	auto             local_data() const&{return local_cutout().data_elements();}
	multi::extensions_type_<2> extensions() const&{return {n0_, std::get<1>(local_cutout().extensions())};}
	friend auto extensions(array const& self){return self.extensions();}

	array& operator=(multi::array<T, 2> const& other) &{
		if(other.extensions() == extensions()) local_cutout() = other.stenciled(std::get<0>(local_cutout().extensions()), std::get<1>(local_cutout().extensions()));
		else{
			array tmp{other};
			std::swap(*this, tmp);
		}
		return *this;
	}
	template<class Array, class=std::enable_if_t<not std::is_same<Array, multi::array<T, 2>>{}> >
	array& operator=(Array const& other) &{
		assert( other.extensions() == this->extensions() );
		
		static_assert( std::is_trivially_assignable<T&, T>{}, "!" );
		static_assert( sizeof(T)%sizeof(double)==0, "!" );
		
		auto options = FFTW_ESTIMATE;
		if(other.layout_.is_transposed){
			options |= FFTW_MPI_TRANSPOSED_IN;
			n0_ = std::get<1>(other.extensions()).size();
		}
		
		fftw_plan p = fftw_mpi_plan_many_transpose(
			std::get<0>(extensions()).size(), std::get<1>(extensions()).size(), sizeof(T)/sizeof(double), 
			FFTW_MPI_DEFAULT_BLOCK, other.layout_.block,
			reinterpret_cast<double*>(const_cast<T*>(other.local_cutout().base())), 
			reinterpret_cast<double*>(this->local_cutout().data_elements()),
			this->comm_.get(), options
		);
		fftw_execute(p);
		fftw_destroy_plan(p);
		
		local_ptr_ = array_ptr<T, 2, local_pointer_t>{this->local_cutout().data_elements(), local_extension_2d(other.extensions())};
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
	basic_array<T, typename std::pointer_traits<local_pointer_t>::template rebind<T const>> transposed() const{
		return basic_array<T, typename std::pointer_traits<local_pointer_t>::template rebind<T const>>{
			layout_t{n0_, true, FFTW_MPI_DEFAULT_BLOCK}, this->local_cutout().layout().transpose(), this->local_cutout().data_elements()
		};
	}
};

#endif
#endif

#if 0
boost::multi::fftw::mpi::scattered_array<std::complex<double>, 2>& dft(
	boost::multi::fftw::mpi::scattered_array<std::complex<double>, 2> const& A, 
	boost::multi::fftw::mpi::scattered_array<std::complex<double>, 2>      & B, 
	fftw::sign /*s*/
){
	(void)A;
//	assert( A.extensions() == B.extensions() );
//	assert( A.comm() == B.comm() );
#if 0
	fftw_plan p = fftw_mpi_plan_dft_2d(
		std::get<0>(A.extensions()).size(), std::get<1>(A.extensions()).size(), 
		(fftw_complex *)A.local_cutout().data_elements(), (fftw_complex *)B.local_cutout().data_elements(), 
		A.comm().get(),
		s, FFTW_ESTIMATE
	);
	fftw_execute(p);
	fftw_destroy_plan(p);
#endif
	return B;
}
#endif

#if 0
array_transposed<std::complex<double>, 2>& dft(
	array<std::complex<double>, 2> const& A, 
	array_transposed<std::complex<double>, 2>& B, 
	fftw::sign s
){
// http://www.fftw.org/fftw3_doc/MPI-Plan-Creation.html
//	assert( A.extensions() == B.extensions() );
	assert( A.comm() == B.comm() );
	fftw_plan p = fftw_mpi_plan_dft_2d(
		std::get<0>(A.extensions()).size(), std::get<1>(A.extensions()).size(), 
		(fftw_complex *)A.local_cutout().data_elements(), (fftw_complex *)B.local_cutout().data_elements(), 
		A.comm().get(),
		s, FFTW_ESTIMATE | FFTW_MPI_TRANSPOSED_OUT
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
#endif

#if not __INCLUDE_LEVEL__

#include<boost/mpi3/main_environment.hpp>
#include<boost/mpi3/ostream.hpp>


#include "../../fftw/mpi/environment.hpp"

namespace mpi3 = boost::mpi3;
namespace multi = boost::multi;
namespace fftw = multi::fftw;
namespace mpi = fftw::mpi;

int mpi3::main(int, char*[], mpi3::environment& menv){
	mpi::environment fenv(menv);

	auto world = menv.world();
	mpi3::ostream os{world};
	
	using T = std::complex<double>;

	mpi::scattered_array<T, 2> S({14, 19}, world);
	
	using std::get;
	if(auto x = extensions(S.local_cutout()))
		for(auto i : get<0>(x))
			for(auto j : get<1>(x))
				S.local_cutout()[i][j] = T(i, j);//std::complex<double>(i + j, i + 2*j + 1)/std::abs(std::complex<double>(i + j, i + 2*j + 1));

	mpi::gathered_array<T, 2> G = S.gather();
	G.local_cutout();

	assert( G.extensions() == {14, 19} );
	if(world.rank() == 0){
		assert( G.extensions() == {14, 19} );
		assert( G.local_cutout().extensions() == {14, 19} );
	}
	if(world.rank() != 0){
		assert( G.extensions() == {14, 19} );
		assert( G.local_cutout().extensions() == {0, 0} );
	}

	multi::array<T, 2> A = S.gather();
	if(world.rank() == 0) assert( A.extensions() == {14, 19} );
	if(world.rank() != 0) assert( A.empty() );
		
	world.barrier();
	if(world.root()){
		std::cout<<"-------------\n";
		if(auto x = extensions(G.local_cutout()))
			for(auto i : get<0>(x)){
				for(auto j : get<1>(x))
					std::cout<< G.local_cutout()[i][j] <<'\t';
				std::cout<<std::endl;
			}
	}else assert(G.local_cutout().empty());
	
	mpi::scattered_array<T, 2> S2 = G.scatter();

	assert( S2 == S );


	mpi::gathered_array<T, 2> G2 = S2.gather();
	
	if(world.root()){
		std::cout<<"-------------\n";
		if(auto x = extensions(G2.local_cutout()))
			for(auto i : get<0>(x)){
				for(auto j : get<1>(x))
					std::cout<< G2.local_cutout()[i][j] <<'\t';
				std::cout<<std::endl;
			}
		assert( G2.local_cutout() == G.local_cutout() );
	}else assert(G2.local_cutout().empty());

//	assert( S == S2 );

//	if(not world.root()) assert( G.local_cutout().empty() );
	
//	mpi::gathered_array<double, 2> G({8, 15}, world);

/*	
	auto const A = [&]{
		os<<"global sizes"<< std::get<0>(A.extensions()) <<'x'<< std::get<1>(A.extensions()) <<' '<< A.num_elements() <<std::endl;
		os<< A.local_cutout().extension() <<'x'<< std::get<1>(A.local_cutout().extensions()) <<"\t#="<< A.local_cutout().num_elements() <<" allocated "<< A.local_count() <<std::endl;
		if(auto x = A.local_cutout().extensions())
			for(auto i : std::get<0>(x))
				for(auto j : std::get<1>(x))
					A.local_cutout()[i][j] = i + j;//std::complex<double>(i + j, i + 2*j + 1)/std::abs(std::complex<double>(i + j, i + 2*j + 1));
		return A;
	}();
*/
/*
	multi::fftw::mpi::scattered_array<std::complex<double>, 2> B(A.extensions(), world);
	
	multi::array<std::complex<double>, 2> A2 = A;
	assert( A2 == A );
	
	using multi::fftw::dft_forward;
*/
#if 0
	dft_forward(A , B );
	dft_forward(A2, A2);

	{
		auto x = B.local_cutout().extensions();
		for(auto i : std::get<0>(x))
			for(auto j : std::get<1>(x))
				if(not( std::abs(B.local_cutout()[i][j] - A2[i][j]) < 1e-12 )){
					std::cout<< B.local_cutout()[i][j] - A2[i][j] <<' '<< std::abs(B.local_cutout()[i][j] - A2[i][j]) <<'\n';
				}
	}
	
	multi::fftw::mpi::array_transposed<std::complex<double>, 2> AT(A.extensions(), world);
	os<< "global sizes" << std::get<0>(AT.extensions()) <<'x'<< std::get<1>(AT.extensions()) <<' '<< AT.num_elements() <<std::endl;
	os<< AT.local_cutout().extension() <<'x'<< std::get<1>(AT.local_cutout().extensions()) <<"\t#="<< AT.local_cutout().num_elements() <<" allocated "<< AT.local_count() <<std::endl;

	dft(A, AT, multi::fftw::forward);
	
	if(world.rank() == 0){
		if(auto x = B.local_cutout().extensions()){
			for(auto i : std::get<0>(x)){
				for(auto j : std::get<1>(x))
					std::cout<< B.local_cutout()[i][j] <<' ';
				std::cout<<'\n';
			}
		}
		
		if(auto x = AT.local_cutout().extensions()){
			for(auto i : std::get<0>(x)){
				for(auto j : std::get<1>(x))
					std::cout<< AT.local_cutout()[i][j] <<' ';
				std::cout<<'\n';
			}
		}
	}
#endif
	return 0;
}
#endif
#endif

