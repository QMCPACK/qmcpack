// Copyright 2021-2026 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

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

#endif
