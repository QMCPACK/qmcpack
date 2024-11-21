// Copyright 2019-2023 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_MULTI_ADAPTORS_BLAS_TRSM_HPP
#define BOOST_MULTI_ADAPTORS_BLAS_TRSM_HPP

#include <boost/multi/adaptors/blas/core.hpp>
#include <boost/multi/adaptors/blas/filling.hpp>
#include <boost/multi/adaptors/blas/operations.hpp>
#include <boost/multi/adaptors/blas/side.hpp>

namespace boost::multi::blas {

enum class diagonal : char {
	    unit = 'U',
	non_unit = 'N', general = non_unit
};

template<blas::filling Fill, class Array>
class triangular_part {
	Array const& ref_;  // NOLINT(cppcoreguidelines-avoid-const-or-ref-data-members)

	public:
	explicit triangular_part(Array const& ref) : ref_{ref} {}
	static constexpr auto filling() { return Fill; }
	using underlying_type = Array;
	auto underlying() const -> Array const& { return ref_;}
};

template<blas::filling Fill, class Array>
auto triangular_parted(Array const& arr) {
	return triangular_part<Fill, Array>{arr};
}

template<class Array>
auto lower_parted(Array const& arr) {return triangular_parted<filling::lower>(arr);}

template<class Array>
auto upper_parted(Array const& arr) {return triangular_parted<filling::upper>(arr);}

template<class Array> auto L(Array const& arr) { return lower_parted(arr); }  // NOLINT(readability-identifier-naming) BLAS naming
template<class Array> auto U(Array const& arr) { return upper_parted(arr); }  // NOLINT(readability-identifier-naming) BLAS naming

template<class Matrix>
auto triangular(multi::blas::filling f, Matrix const& m) {  // NOLINT(readability-identifier-length) BLAS naming
	auto ret =+ m;
	switch(f) {
	case multi::blas::filling::upper:
		{
			auto ext = extension(ret);
			std::for_each(ext.begin(), ext.end(), [&ret](auto idx) {
				std::fill_n(ret[idx].begin(), std::min(idx, size(~ret)), 0.0);
			});
		}
		break;
	case multi::blas::filling::lower:
		{
			auto extt = extension(~ret);
			std::for_each(extt.begin(), extt.end(), [&ret](auto jdx) {
				std::fill_n( (~ret)[jdx].begin(), std::min(jdx, size( ret)), 0.0);
			});
		}
		break;
	}
	return ret;
}

using core::trsm;

template<class Context, class A2D, class B2D>
auto trsm(Context&& ctxt, blas::side a_side, blas::filling a_fill, blas::diagonal a_diag, typename A2D::element_type alpha, A2D const& a, B2D&& b) // NOLINT(readability-function-cognitive-complexity,readability-identifier-length) cognitive load 115, BLAS naming
-> B2D&& {
	if(a_side == blas::side::left ) {assert(size(~a) >= size( b));}  // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay)
	if(a_side == blas::side::right) {assert(size( a) >= size(~b));}  // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay)

	assert( stride( a) == 1 || stride(~a) == 1 );  // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay)
	assert( stride( b) == 1 || stride(~b) == 1 );  // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay)

	if(size(b) != 0) {
		#define CTXT std::forward<Context>(ctxt)
		if       constexpr(! is_conjugated<A2D>{} && ! is_conjugated<B2D>{}) {
			if     (stride( a)==1 && stride( b)==1) {CTXT->trsm(static_cast<char>(    (a_side)), static_cast<char>(-a_fill), 'N', static_cast<char>(a_diag), size( b), size(~b),      alpha ,            a.base() , stride(~a),            b.base() , stride(~b));}
			else if(stride(~a)==1 && stride(~b)==1) {CTXT->trsm(static_cast<char>(swap(a_side)), static_cast<char>(+a_fill), 'N', static_cast<char>(a_diag), size(~b), size( b),      alpha ,            a.base() , stride( a),            b.base() , stride( b));}
			else if(stride( a)==1 && stride(~b)==1) {CTXT->trsm(static_cast<char>(swap(a_side)), static_cast<char>(-a_fill), 'T', static_cast<char>(a_diag), size(~b), size( b),      alpha ,            a.base() , stride(~a),            b.base() , stride( b));}
			else if(stride(~a)==1 && stride( b)==1) {CTXT->trsm(static_cast<char>(    (a_side)), static_cast<char>(+a_fill), 'T', static_cast<char>(a_diag), size( b), size(~b),      alpha ,            a.base() , stride( a),            b.base() , stride(~b));}
			else                                    {assert(0 && "not implemented in blas");}  // LCOV_EXCL_LINE
		} else if constexpr(   is_conjugated<A2D>{} && ! is_conjugated<B2D>{}) {
			if     (stride( a)==1 && stride(~b)==1) {CTXT->trsm(static_cast<char>(swap(a_side)), static_cast<char>(-a_fill), 'C', static_cast<char>(a_diag), size(~b), size( b),      alpha , underlying(a.base()), stride(~a),            b.base() , stride( b));}
			else if(stride(~a)==1 && stride( b)==1) {CTXT->trsm(static_cast<char>(    (a_side)), static_cast<char>(+a_fill), 'C', static_cast<char>(a_diag), size( b), size(~b),      alpha , underlying(a.base()), stride( a),            b.base() , stride(~b));}
			else                                    {assert(0 && "not implemented in blas");}  // LCOV_EXCL_LINE
		} else if constexpr(! is_conjugated<A2D>{} &&    is_conjugated<B2D>{}) {
			if     (stride(~a)==1 && stride( b)==1) {CTXT->trsm(static_cast<char>(    (a_side)), static_cast<char>(+a_fill), 'C', static_cast<char>(a_diag), size( b), size(~b), conj(alpha),            a.base() , stride( a), underlying(b.base()), stride(~b));}
		//  else if(stride( a)==1 && stride(~b)==1) {assert(0 && "not implemented in blas");}  // LCOV_EXCL_LINE
			else                                    {assert(0 && "not implemented in blas");}  // LCOV_EXCL_LINE
		} else if constexpr(   is_conjugated<A2D>{} &&     is_conjugated<B2D>{}) {
			if     (stride( a)==1 && stride(~b)==1) {CTXT->trsm(static_cast<char>(swap(a_side)), static_cast<char>(-a_fill), 'T', static_cast<char>(a_diag), size(~b), size( b), conj(alpha), underlying(a.base()), stride(~a), underlying(b.base()), stride( b));}
			else if(stride(~a)==1 && stride( b)==1) {CTXT->trsm(static_cast<char>(    (a_side)), static_cast<char>(+a_fill), 'T', static_cast<char>(a_diag), size( b), size(~b), conj(alpha), underlying(a.base()), stride( a), underlying(bbase(b)), stride(~b));}
			else                                    {assert(0 && "not implemented in blas");}  // LCOV_EXCL_LINE
		}
		#undef CTXT
	}
	return std::forward<B2D>(b);
}

template<class A2D, class B2D>
auto trsm(blas::side a_side, blas::filling a_fill, blas::diagonal a_diag, typename A2D::element_type alpha, A2D const& a, B2D&& b) -> decltype(auto) {  // NOLINT(readability-identifier-length) BLAS naming
	if constexpr(! is_conjugated<A2D>{}) {return trsm(default_context_of(           a.base() ), a_side, a_fill, a_diag, alpha, a, std::forward<B2D>(b));}
	else                                 {return trsm(default_context_of(underlying(a.base())), a_side, a_fill, a_diag, alpha, a, std::forward<B2D>(b));}
}

template<class Context, class A2D, class B2D>
auto trsm(Context&& ctxt, blas::side a_side, blas::filling a_fill, typename A2D::element_type alpha, A2D const& a, B2D&& b)  // NOLINT(readability-identifier-length) BLAS naming
->decltype(trsm(std::forward<Context>(ctxt), a_side, a_fill, blas::diagonal::non_unit, alpha, a, std::forward<B2D>(b))) {
	return trsm(std::forward<Context>(ctxt), a_side, a_fill, blas::diagonal::non_unit, alpha, a, std::forward<B2D>(b)); }

#if defined __NVCC__
	#ifdef __NVCC_DIAG_PRAGMA_SUPPORT__
		#pragma nv_diagnostic push
		#pragma nv_diag_suppress = implicit_return_from_non_void_function
	#else
		#pragma    diagnostic push
		#pragma    diag_suppress = implicit_return_from_non_void_function
	#endif
#elif defined __NVCOMPILER
	#pragma    diagnostic push
	#pragma    diag_suppress = implicit_return_from_non_void_function
#endif
template<class A2D, class B2D>
auto trsm(blas::side a_side, blas::filling a_fill, typename A2D::element_type alpha, A2D const& a, B2D&& b) -> decltype(auto) {  // NOLINT(readability-identifier-length) BLAS naming
	if constexpr(! is_conjugated<A2D>{}) {return trsm(blas::default_context_of(           a.base() ), a_side, a_fill, alpha, a, std::forward<B2D>(b));}
	else                                 {return trsm(blas::default_context_of(underlying(a.base())), a_side, a_fill, alpha, a, std::forward<B2D>(b));}
}
#if defined __NVCC__
	#ifdef __NVCC_DIAG_PRAGMA_SUPPORT__
		#pragma nv_diagnostic pop
	#else
		#pragma    diagnostic pop
	#endif
#elif defined __NVCOMPILER
	#pragma    diagnostic pop
#endif

template<class UTArr, class B2D>
auto trsm(blas::side a_side, typename UTArr::underlying_type::element_type alpha, UTArr const& a, B2D&& b)  // NOLINT(readability-identifier-length) BLAS naming
->decltype(trsm(a_side, a.filling(), blas::diagonal::non_unit, alpha, a.underlying(), std::forward<B2D>(b))) {
	return trsm(a_side, a.filling(), blas::diagonal::non_unit, alpha, a.underlying(), std::forward<B2D>(b)); }

namespace operators {

	template<class B2D, class UL>
	auto operator/=(B2D&& b, UL const& a)  // NOLINT(readability-identifier-length) BLAS naming
	->decltype(blas::trsm(blas::side::right, 1.0, a, std::forward<B2D>(b))) {
		return blas::trsm(blas::side::right, 1.0, a, std::forward<B2D>(b)); }

	template<class B2D, class UL>
	auto operator|=(B2D&& b, UL const& a)  // NOLINT(readability-identifier-length) BLAS naming
	->decltype(blas::trsm(blas::side::left, 1.0, a, std::forward<B2D>(b))) {
		return blas::trsm(blas::side::left, 1.0, a, std::forward<B2D>(b)); }

	using blas::U;
	using blas::L;

}  // end namespace operators

}  // end namespace boost::multi::blas

#endif
