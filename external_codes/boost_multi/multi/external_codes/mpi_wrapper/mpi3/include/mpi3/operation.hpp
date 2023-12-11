// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2018-2023 Alfredo A. Correa

#ifndef BOOST_MPI3_OPERATION_HPP
#define BOOST_MPI3_OPERATION_HPP

#include <mpi3/detail/datatype.hpp>
#include <mpi3/handle.hpp>


#if defined(__NVCC__)
#include <thrust/transform.h>
#endif

#include <algorithm>  // std::transform_n
#include <utility> // std::forward

namespace boost {
namespace mpi3 {

template<class T>
struct commutative_operation;

template<class T>
struct builtin_operation;

template<class T, class P = T*>
struct operation : detail::nondefault_handle<operation<T>, MPI_Op, MPI_Op_free> { // user_operation, operator_ . operator is a C++ keyword
	using base = typename detail::nondefault_handle<operation<T>, MPI_Op, MPI_Op_free>;
	using typename base::nondefault_handle;
	operation() = delete;

	template<class Op>
	// NOLINTNEXTLINE(bugprone-easily-swappable-parameters) signature is fixed
	constexpr static void combine(void const* in, void* inout, int const* len, MPI_Datatype* /*dtype*/) {  // cppcheck-suppress constParameter ; signature is fixed
		auto in_t    = reinterpret_cast<typename std::pointer_traits<P>::template rebind<T const>&>(in   );  // NOLINT(llvm-qualified-auto,readability-qualified-auto,cppcoreguidelines-pro-type-reinterpret-cast)
		auto inout_t = reinterpret_cast<                             P                           &>(inout);  // NOLINT(llvm-qualified-auto,readability-qualified-auto,cppcoreguidelines-pro-type-reinterpret-cast)
		#if defined(__NVCC__)
		thrust::transform(
		#else
		std::transform(
		#endif
			in_t, std::next(in_t, *len), inout_t, inout_t,
			[]
			#if defined(__NVCC__)
			__host__ __device__
			#endif
			(T const& a, T const& b)
			#if not defined(__NVCC__)
			constexpr
			#endif
			{return Op{}(a, b);}
		);
	//	for(int i = 0; i != *len; i++) {
	//		inout_t[i] = Op{}(std::move(inout_t[i]), in_t[i]);  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
	//	}
	}

	public:
	template<class F>
	explicit operation(F const& /*f*/) : base(detail::uninitialized{}) {
		MPI_Op_create(reinterpret_cast<MPI_User_function*>(&combine<F>), /*commutative*/ true, &(this->impl_));  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast)
	}

	template<class F, typename = std::enable_if_t<not std::is_same<std::decay_t<F>, operation>{}> >
	operation(F&& f, bool commutative) : base(detail::uninitialized{}) {
		MPI_Op_create(
			&f,
		//	reinterpret_cast<void (*)(void*, void*, int*, int*)>(&f),
			commutative,
			&(this->impl_)
		);
	}

	operation(operation const&) = delete;
	operation(operation     &&) = delete;

	operation& operator=(operation const&) = delete;
	operation& operator=(operation     &&) = delete;

	~operation() = default;

	MPI_Op operator&() const {return this->impl_;}  // NOLINT(google-runtime-operator)

#if 0
	enum struct code : MPI_Op{
		maximum = MPI_MAX, minimum = MPI_MIN,
		sum = MPI_SUM, product = MPI_PROD,
		logical_and = MPI_LAND, bitwise_and = MPI_BAND,
		logical_or  = MPI_LOR,   bitwise_or = MPI_BOR,
		logical_xor = MPI_LXOR, bitwise_xor = MPI_BXOR,
		max_value_location = MPI_MAXLOC,
		min_value_location = MPI_MINLOC
	};
#endif
//	operation(operation::code c) : base((MPI_Op)c){}

//	static operation const sum;//(operation::code::sum);
//	static operation const product;
//	static operation const maximum;
//	static operation const minimum;

};

//operation const sum    (operation::code::sum);
//operation const product(operation::code::product);
//operation const maximum(operation::code::maximum);
//operation const minimum(operation::code::minimum);

template<class T = void>
using plus = std::plus<T>;
template<class T = void>
using minus = std::minus<T>;
template<class T = void>
using multiplies = std::multiplies<T>;

template<class T = void> struct min {
	constexpr T const& operator()(T const& t1, T const& t2) const {return std::min(t1, t2);}
};
template<> struct min<void>{
	template<class T1, class T2> constexpr decltype(auto) operator()(T1&& t1, T2&& t2) const {return std::min(std::forward<T1>(t1), std::forward<T2>(t2));}
};

template<class T = void> struct max {
	constexpr T const& operator()(T const& t1, T const& t2) const {return std::max(t1, t2);}
};
template<> struct max<void> {
	template<class T1, class T2> constexpr decltype(auto) operator()(T1&& t1, T2&& t2) const {return std::max(std::forward<T1>(t1), std::forward<T2>(t2));}
};

template<class T = void> struct max_loc {  // the only differences is that argument is assumed to be a pair, and second element is int
	constexpr T const& operator()(T const& t1, T const& t2) const {std::max(t1, t2);}
};
template<> struct max_loc<void> {
	template<class T1, class T2> constexpr decltype(auto) operator()(T1&& t1, T2&& t2) const {return std::max(std::forward<T1>(t1), std::forward<T2>(t2));}
};

template<class T = void> struct min_loc {  // the only differences is that argument is assumed to be a pair, and second element is int
	constexpr T const& operator()(T const& t1, T const& t2) const {std::min(t1, t2);}
};
template<> struct min_loc<void> {
	template<class T1, class T2> constexpr decltype(auto) operator()(T1&& t1, T2&& t2) const {return std::min(std::forward<T1>(t1), std::forward<T2>(t2));}
};

template<>
struct operation<double> {
 private:
	MPI_Op impl_;

 public:
	explicit operation(std::plus<> /*op*/) : impl_{MPI_SUM} {}
	explicit operation(std::plus<double> /*op*/) : impl_{MPI_SUM} {}

	explicit operation(mpi3::min<> /*op*/) : impl_{MPI_MIN} {}
	explicit operation(mpi3::min<double> /*op*/) : impl_{MPI_MIN} {}

	explicit operation(mpi3::max<> /*op*/) : impl_{MPI_MAX} {}
	explicit operation(mpi3::max<double> /*op*/) : impl_{MPI_MAX} {}

	MPI_Op operator&() const {return impl_;}  // NOLINT(google-runtime-operator)
};

template<>
struct operation<mpi3::vlp<int> > {
 private:
	MPI_Op impl_;

 public:
	// explicit operation(std::plus<> /*op*/) : impl_{MPI_SUM} {}
	// explicit operation(std::plus<double> /*op*/) : impl_{MPI_SUM} {}

	// explicit operation(mpi3::min<> /*op*/) : impl_{MPI_MIN} {}
	// explicit operation(mpi3::min<double> /*op*/) : impl_{MPI_MIN} {}

	explicit operation(mpi3::max_loc<> /*op*/) : impl_{MPI_MAXLOC} {}
	explicit operation(mpi3::max_loc<mpi3::vlp<int> > /*op*/) : impl_{MPI_MAXLOC} {}

	MPI_Op operator&() const {return impl_;}  // NOLINT(google-runtime-operator)
};


template<class Op> struct predefined_operation;

// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define BOOST_MPI3_DECLARE_PREDEFINED_OPERATION(CppoP, MpinamE, NamE) \
template<> struct predefined_operation<CppoP>{ \
/*	constexpr*/ operator MPI_Op() const{return MpinamE;} \
/*	static constexpr MPI_Op value = MpinamE;*/ \
}; \
using NamE = predefined_operation<CppoP>  // NOLINT(bugprone-macro-parentheses)

BOOST_MPI3_DECLARE_PREDEFINED_OPERATION(std::plus<>       , MPI_SUM , sum        );
BOOST_MPI3_DECLARE_PREDEFINED_OPERATION(std::multiplies<> , MPI_PROD, product    );
BOOST_MPI3_DECLARE_PREDEFINED_OPERATION(std::logical_and<>, MPI_LAND, logical_and);

BOOST_MPI3_DECLARE_PREDEFINED_OPERATION(std::bit_and<>, MPI_BAND, bitwise_and);
BOOST_MPI3_DECLARE_PREDEFINED_OPERATION(std::bit_or<> , MPI_BOR , bitwise_or);

BOOST_MPI3_DECLARE_PREDEFINED_OPERATION(max<>, MPI_MAX, maximum);
BOOST_MPI3_DECLARE_PREDEFINED_OPERATION(min<>, MPI_MIN, minimum);

BOOST_MPI3_DECLARE_PREDEFINED_OPERATION(max_loc<>, MPI_MAXLOC, maximum_location);
BOOST_MPI3_DECLARE_PREDEFINED_OPERATION(min_loc<>, MPI_MINLOC, minimum_location);

#undef BOOST_MPI3_DECLARE_PREDEFINED_OPERATION

template<class T>
struct commutative_operation : operation<T> {
	template<class F,  typename = std::enable_if_t<not std::is_same<std::decay_t<F>, operation<T> >{}> >
	explicit commutative_operation(F&& f) : operation<T>(std::forward<F>(f), true){}
};

template<class T>
struct non_commutative_operation : operation<T> {
	template<class F,  typename = std::enable_if_t<not std::is_same<std::decay_t<F>, operation<T> >{}>>
	explicit non_commutative_operation(F&& f) : operation<T>(std::forward<F>(f), false){}
};

}  // end namespace mpi3
}  // end namespace boost

#endif
