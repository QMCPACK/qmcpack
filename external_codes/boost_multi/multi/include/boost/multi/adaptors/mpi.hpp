// Copyright 2024-2026 Alfredo A. Correa
// Distributed under the Boost Software License, Version 10.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_MULTI_ADAPTORS_MPI_HPP
#define BOOST_MULTI_ADAPTORS_MPI_HPP

#include "boost/multi/array.hpp"

#include <type_traits>

#define OMPI_SKIP_MPICXX 1
#include <mpi.h>

#include <cassert>  // for assert
#include <complex>
#include <limits>   // for numeric_limits  NOLINT(misc-include-cleaner)
#include <utility>  // for exchange, move

namespace boost::multi::mpi {

using const_MPI_Datatype = MPI_Datatype const;

template<class T> static inline const_MPI_Datatype const datatype = std::conditional_t<static_cast<bool>(sizeof(T*)), void, int>{};

// template<> MPI_Datatype const datatype<char> = MPI_CHAR;
// template<> MPI_Datatype const datatype<unsigned char> = MPI_UNSIGNED_CHAR;

// #if(__cplusplus >= 201703L)
// MPI3_DECLARE_DATATYPE(byte                   , MPI_BYTE);
// #endif
// MPI3_DECLARE_DATATYPE(wchar                  , MPI_WCHAR);

// MPI3_DECLARE_DATATYPE(short                  , MPI_SHORT);
// MPI3_DECLARE_DATATYPE(unsigned short         , MPI_UNSIGNED_SHORT);
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wmissing-variable-declarations"
#endif

template<> inline const_MPI_Datatype const datatype<int> = MPI_INT;  // NOLINT(misc-misplaced-const)

// MPI3_DECLARE_DATATYPE(unsigned int           , MPI_UNSIGNED);
// MPI3_DECLARE_DATATYPE(long                   , MPI_LONG);
// MPI3_DECLARE_DATATYPE(unsigned long          , MPI_UNSIGNED_LONG);
// MPI3_DECLARE_DATATYPE(float                  , MPI_FLOAT);

template<> inline const_MPI_Datatype const datatype<float>  = MPI_FLOAT;

template<> inline const_MPI_Datatype const datatype<double> = MPI_DOUBLE;

template<> inline const_MPI_Datatype const datatype<std::complex<double>> = MPI_DOUBLE_COMPLEX;

#ifdef __clang__
#pragma clang diagnostic pop
#endif

// MPI3_DECLARE_DATATYPE(long double            , MPI_LONG_DOUBLE);
// MPI3_DECLARE_DATATYPE(long long int          , MPI_LONG_LONG_INT);

// MPI3_DECLARE_DATATYPE(bool                   , MPI_C_BOOL);  // C++ binding not used MPI_CXX_BOOL);

template<typename T>
struct datatype_t {
	static const_MPI_Datatype const value;  // = datatype<T>;
	auto operator()() const -> decltype(datatype<T>) { return datatype<T>; }
	operator const_MPI_Datatype() const { return datatype<T>; }  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions,readability-const-return-type)
};

template<typename T>
const_MPI_Datatype const datatype_t<T>::value = datatype<T>;

// class data {
// 	void*        buf_;
// 	MPI_Datatype datatype_;

//  public:
// 	data(void* buf, MPI_Datatype datatype) : buf_{buf}, datatype_{datatype} {}
// 	template<class It>
// 	explicit data(It first)                                            // NOLINT(cppcoreguidelines-pro-type-member-init,hicpp-member-init)
// 	: buf_{const_cast<void*>(static_cast<void const*>(first.base()))}  // NOLINT(cppcoreguidelines-pro-type-const-cast)
// 	{
// 		assert( first.stride() <= std::numeric_limits<int>::max() );
// 		MPI_Type_vector(
// 			1, 1,
// 			static_cast<int>(first.stride()),
// 			mpi::datatype<typename It::element>,
// 			&datatype_
// 		);

// 		MPI_Type_commit(&datatype_);  // type cannot be used until committed, in communication operations at least
// 	}

// 	data(data const&) = delete;
// 	data(data&&)      = delete;

// 	auto operator=(data const&) = delete;
// 	auto operator=(data&&)      = delete;

// 	~data() { MPI_Type_free(&datatype_); }

// 	auto buffer() const { return buf_; }
// 	auto datatype() const { return datatype_; }
// };

template<class Layout>
auto create_subarray_aux(
	Layout        lyt,
	int           subcount,
	MPI_Datatype  old_datatype,
	MPI_Datatype* new_datatype
) -> int {
	MPI_Datatype sub_type;  // NOLINT(cppcoreguidelines-init-variables)

	if constexpr(Layout::dimensionality == 1) {
		MPI_Type_dup(old_datatype, &sub_type);
	} else {
		create_subarray_aux(lyt.sub(), lyt.sub().size(), old_datatype, &sub_type);
	}

	int dt_size;  // NOLINT(cppcoreguidelines-init-variables)
	MPI_Type_size(old_datatype, &dt_size);

	assert(lyt.stride() * dt_size <= std::numeric_limits<MPI_Aint>::max());
	{
		MPI_Datatype vector_datatype;  // NOLINT(cppcoreguidelines-init-variables)
		MPI_Type_create_hvector(
			subcount, 1,
			lyt.stride() * dt_size,
			sub_type, &vector_datatype
		);

		MPI_Type_create_resized(vector_datatype, 0, lyt.stride() * dt_size, new_datatype);
		MPI_Type_free(&vector_datatype);
	}
	MPI_Type_free(&sub_type);
	return MPI_SUCCESS;
}

template<class T = void, template<typename> class DatatypeT = mpi::datatype_t, class Size = int>
class skeleton {
	MPI_Datatype datatype_ = MPI_DATATYPE_NULL;
	Size         count_ = 0;

	skeleton() = default;

	auto operator=(skeleton&& other) & noexcept -> skeleton& {
		count_    = other.count_;
		datatype_ = std::exchange(other.datatype_, MPI_DATATYPE_NULL);
		return *this;
	}

 public:
	template<class Stride, class SubLayout>
	skeleton(Stride stride, SubLayout const& sublyt, MPI_Datatype dt, Size subcount) : count_{1} {
		MPI_Datatype              sub_type;  // NOLINT(cppcoreguidelines-init-variables)
		[[maybe_unused]] skeleton sk;        // NOLINT(misc-const-correctness)
		if constexpr(SubLayout::dimensionality == 0) {
			sub_type = dt;
		} else {
			assert( sublyt.size() <= std::numeric_limits<int>::max() );
			sk       = skeleton(sublyt, dt, static_cast<int>(sublyt.size()));
			sub_type = sk.datatype();
		}

		int dt_size;  // NOLINT(cppcoreguidelines-init-variables)
		MPI_Type_size(dt, &dt_size);

		{
			MPI_Datatype vector_datatype;  // NOLINT(cppcoreguidelines-init-variables)
			MPI_Type_create_hvector(
				subcount, 1,
				stride * dt_size,
				sub_type, &vector_datatype
			);

			MPI_Type_create_resized(vector_datatype, 0, stride * dt_size, &datatype_);
			MPI_Type_free(&vector_datatype);
		}
	}

	template<class Layout>
	skeleton(Layout const& lyt, MPI_Datatype dt, Size subcount) : count_{static_cast<Size>(lyt.size())} {
		assert(lyt.size() <= std::numeric_limits<Size>::max());

		MPI_Datatype              sub_type;  // NOLINT(cppcoreguidelines-init-variables)
		[[maybe_unused]] skeleton sk;        // NOLINT(misc-const-correctness)
		if constexpr(Layout::dimensionality == 1) {
			sub_type = dt;
		} else {
			assert( lyt.sub().size() <= std::numeric_limits<int>::max() );
			sk       = skeleton(lyt.sub(), dt, static_cast<int>(lyt.sub().size()));
			sub_type = sk.datatype();
		}

		int dt_size;  // NOLINT(cppcoreguidelines-init-variables)
		MPI_Type_size(dt, &dt_size);

		{
			MPI_Datatype vector_datatype;  // NOLINT(cppcoreguidelines-init-variables)
			MPI_Type_create_hvector(
				subcount, 1,
				lyt.stride() * dt_size,
				sub_type, &vector_datatype
			);

			MPI_Type_create_resized(vector_datatype, 0, lyt.stride() * dt_size, &datatype_);
			MPI_Type_free(&vector_datatype);
		}
	}

 	skeleton(skeleton&& other) noexcept
	: count_{other.count_}, datatype_{std::exchange(other.datatype_, MPI_DATATYPE_NULL)} {}

	template<class Stride, class SubLayout>
	skeleton(Stride stride, SubLayout const& sublyt, MPI_Datatype dt)
	: skeleton{stride, sublyt, dt, 1} {
		MPI_Type_commit(&datatype_);
	}

	template<class Layout>
	skeleton(Layout const& lyt, MPI_Datatype dt)
	: skeleton{lyt, dt, 1} {
		MPI_Type_commit(&datatype_);
	}

	template<class Layout>
	explicit skeleton(Layout const& lyt) : skeleton{lyt, DatatypeT<T>{}()} {}

	skeleton(skeleton const&) = delete;

	auto operator=(skeleton const&) = delete;

	~skeleton() {
		if(datatype_ != MPI_DATATYPE_NULL) {
			MPI_Type_free(&datatype_);
		}
	}

	auto count() const { return count_; }
	auto datatype() const& { return datatype_; }
	auto datatype() && { return std::exchange(datatype_, MPI_DATATYPE_NULL); }
};

#if defined(__cpp_deduction_guides) && (__cpp_deduction_guides>=201703L)
template<class Layout> skeleton(Layout, decltype(MPI_INT)) -> skeleton<Layout>;
#endif

template<class Layout>
auto create_subarray(Layout const& lyt, MPI_Datatype old_datatype, MPI_Datatype* new_datatype) -> int {
	int old_datatype_size;  // NOLINT(cppcoreguidelines-init-variables)
	MPI_Type_size(old_datatype, &old_datatype_size);

	// return create_subarray_aux(lyt, 1, old_datatype, new_datatype);
	skeleton const sk(lyt, old_datatype);
	// new_datatype = std::move(sk).type();
	{
		MPI_Datatype vector_datatype;  // NOLINT(cppcoreguidelines-init-variables)
		assert( lyt.size() <= std::numeric_limits<int>::max() );

		MPI_Type_create_hvector(
			static_cast<int>(lyt.size()), 1,
			lyt.stride() * old_datatype_size,
			sk.datatype(), &vector_datatype
		);

		MPI_Type_create_resized(vector_datatype, 0, lyt.stride() * old_datatype_size, new_datatype);
		MPI_Type_free(&vector_datatype);
	}
	return MPI_SUCCESS;
}

template<template<typename> class DatatypeT = mpi::datatype_t, typename Size = int>
class message : skeleton<void, DatatypeT, Size> {
	void* buf_;

	using skeleton_type = skeleton<void, DatatypeT, Size>;

 public:
	message(void* buf, skeleton_type&& sk) : skeleton_type{std::move(sk)}, buf_{buf} {}

	template<class Layout>
	message(void* buf, Layout const& lyt, MPI_Datatype dt) : skeleton_type(lyt, dt), buf_{buf} {}

	template<class ArrayElements>
	explicit message(ArrayElements const& arrelems)
	: message{
		const_cast<void*>(static_cast<void const*>(arrelems.base())),  // NOLINT(cppcoreguidelines-pro-type-const-cast)
		arrelems.layout(),
		DatatypeT<typename ArrayElements::element>{}  // value_type>
	} {}

	message(message const& other) = delete;
	message(message&&)            = delete;

	auto operator=(message const&) = delete;
	auto operator=(message&&)      = delete;

	~message() = default;

	auto buffer() const { return buf_; }
	using skeleton_type::count;
	// auto count() const { return this->count_; }
	using skeleton_type::datatype;
	// auto datatype() const { return this->datatype_; }

	// template<std::size_t Index>
	// std::tuple_element_t<Index, skeleton<>> const& get() const& {
	//  if constexpr(Index == 0)
	//      return buf_;
	//  if constexpr(Index == 1)
	//      return this->count_;
	//  if constexpr(Index == 2)
	//      return this->datatype_;
	// }
};

template<template<typename> class DatatypeT = mpi::datatype_t, typename Size = int>
class iterator : skeleton<void, DatatypeT, Size> {
	void* buf_;
	using skeleton_type = skeleton<void, DatatypeT, Size>;

 public:
	iterator(void* buf, skeleton_type&& sk) : skeleton_type{std::move(sk)}, buf_{buf} {}

	template<class Stride, class SubLayout>
	iterator(void* buf, Stride stride, SubLayout const& sublyt, MPI_Datatype dt) : skeleton_type{stride, sublyt, dt}, buf_{buf} {}

	template<class ArrayIterator, std::enable_if_t<ArrayIterator::rank_v != 1, int> =0>  // NOLINT(modernize-use-constraints) for C++20
	explicit iterator(ArrayIterator const& it)
	: iterator{
		const_cast<void*>(static_cast<void const*>(it.base())),  // NOLINT(cppcoreguidelines-pro-type-const-cast)
		it.stride(),
		it->layout(),
		DatatypeT<typename ArrayIterator::element>{}
	} {}

	template<class ArrayIterator, std::enable_if_t<ArrayIterator::rank_v == 1, int> =0>  // NOLINT(modernize-use-constraints) for C++20
	explicit iterator(ArrayIterator const& it)
	: iterator{
		const_cast<void*>(static_cast<void const*>(it.base())),  // NOLINT(cppcoreguidelines-pro-type-const-cast)
		it.stride(),
		multi::layout_t<0>{},
		DatatypeT<typename ArrayIterator::element>{}
	} {}

	iterator(iterator const& other) = delete;
	iterator(iterator&&)            = delete;

	auto operator=(iterator const&) = delete;
	auto operator=(iterator&&)      = delete;

	~iterator() = default;

	auto buffer() const { return buf_; }
	// auto count() const { return 1; }  // an iterator doesn't have a count in most contexts, you have to figure out the count yourselve, typically size()
	using skeleton_type::count;
	using skeleton_type::datatype;
};

template<template<typename> class DatatypeT = mpi::datatype_t, class Array>
auto begin(Array&& arr) {return iterator<DatatypeT>{std::forward<Array>(arr).begin()}; }

#if defined(__cpp_deduction_guides) && (__cpp_deduction_guides>=201703L)
template<class ArrayElements> message(ArrayElements const&) -> message<>;
#endif

}  // namespace boost::multi::mpi

#endif  // BOOST_MULTI_ADAPTORS_MPI_HPP
