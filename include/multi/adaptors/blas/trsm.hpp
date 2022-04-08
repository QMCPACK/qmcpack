// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2019-2022 Alfredo A. Correa
#pragma once

#include "../blas/core.hpp"
#include "../blas/filling.hpp"
#include "../blas/operations.hpp" // uplo
#include "../blas/side.hpp"

namespace boost::multi::blas {

enum class diagonal : char {
	    unit = 'U', 
	non_unit = 'N', general = non_unit
};

using core::trsm;

template<class Context, class A2D, class B2D>
auto trsm(Context&& ctxt, blas::side a_side, blas::filling a_fill, blas::diagonal a_diag, typename A2D::element_type alpha, A2D const& a, B2D&& b) // NOLINT(readability-function-cognitive-complexity) : cognitive load 115
-> B2D&& {
	if(a_side == blas::side::left ) {assert(size(~a) >= size( b));}  // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay)
	if(a_side == blas::side::right) {assert(size( a) >= size(~b));}  // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay)

	assert( stride( a) == 1 or stride(~a) == 1 );  // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay)
	assert( stride( b) == 1 or stride(~b) == 1 );  // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay)

	if(size(b)!=0){
		#define CTXT std::forward<Context>(ctxt)
		if       constexpr(not is_conjugated<A2D>{} and not is_conjugated<B2D>{}){
			if     (stride( a)==1 and stride( b)==1) {CTXT->trsm(char{    (a_side)}, static_cast<char>(-a_fill), 'N', static_cast<char>(a_diag), size( b), size(~b),      alpha ,            base(a) , stride(~a),            base(b) , stride(~b));}
			else if(stride(~a)==1 and stride(~b)==1) {CTXT->trsm(char{swap(a_side)}, static_cast<char>(+a_fill), 'N', static_cast<char>(a_diag), size(~b), size( b),      alpha ,            base(a) , stride( a),            base(b) , stride( b));}
			else if(stride( a)==1 and stride(~b)==1) {CTXT->trsm(char{swap(a_side)}, static_cast<char>(-a_fill), 'T', static_cast<char>(a_diag), size(~b), size( b),      alpha ,            base(a) , stride(~a),            base(b) , stride( b));}
			else if(stride(~a)==1 and stride( b)==1) {CTXT->trsm(char{    (a_side)}, static_cast<char>(+a_fill), 'T', static_cast<char>(a_diag), size( b), size(~b),      alpha ,            base(a) , stride( a),            base(b) , stride(~b));}
			else                                     {assert(0 && "not implemented in blas");} // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay)
		}else if constexpr(   is_conjugated<A2D>{} and not is_conjugated<B2D>{}){
			if     (stride( a)==1 and stride(~b)==1) {CTXT->trsm(char{swap(a_side)}, static_cast<char>(-a_fill), 'C', static_cast<char>(a_diag), size(~b), size( b),      alpha , underlying(base(a)), stride(~a),            base(b) , stride( b));}
			else if(stride(~a)==1 and stride( b)==1) {CTXT->trsm(char{    (a_side)}, static_cast<char>(+a_fill), 'C', static_cast<char>(a_diag), size( b), size(~b),      alpha , underlying(base(a)), stride( a),            base(b) , stride(~b));}
		//  else if(stride( a)==1 and stride( b)==1) {assert(0 && "not implemented in blas");}  // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay)
		//  else if(stride(~a)==1 and stride(~b)==1) {assert(0 && "not implemented in blas");}  // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay)
			else                                     {assert(0 && "not implemented in blas");}  // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay)
		}else if constexpr(not is_conjugated<A2D>{} and    is_conjugated<B2D>{}){
			if     (stride(~a)==1 and stride( b)==1) {CTXT->trsm(char{    (a_side)}, static_cast<char>(+a_fill), 'C', static_cast<char>(a_diag), size( b), size(~b), conj(alpha),            base(a) , stride( a), underlying(base(b)), stride(~b));}
			else if(stride( a)==1 and stride(~b)==1) {CTXT->trsm(char{swap(a_side)}, static_cast<char>(-a_fill), 'C', static_cast<char>(a_diag), size(~b), size( b), conj(alpha),            base(a) , stride(~a), underlying(base(b)), stride( b));}
		//  else if(stride(~a)==1 and stride(~b)==1) {assert(0 && "not implemented in blas");}  // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay)
		//  else if(stride( a)==1 and stride( b)==1) {assert(0 && "not implemented in blas");}  // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay)
			else                                     {assert(0 && "not implemented in blas");}  // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay)
		}else if constexpr(   is_conjugated<A2D>{} and     is_conjugated<B2D>{}){
			if     (stride( a)==1 and stride(~b)==1) {CTXT->trsm(char{swap(a_side)}, static_cast<char>(-a_fill), 'T', static_cast<char>(a_diag), size(~b), size( b), conj(alpha), underlying(base(a)), stride(~a), underlying(base(b)), stride( b));}
			else if(stride(~a)==1 and stride( b)==1) {CTXT->trsm(char{    (a_side)}, static_cast<char>(+a_fill), 'T', static_cast<char>(a_diag), size( b), size(~b), conj(alpha), underlying(base(a)), stride( a), underlying(base(b)), stride(~b));}
		//  else if(stride(~a)==1 and stride(~b)==1) {assert(0 && "not implemented in blas");}
		//  else if(stride( a)==1 and stride( b)==1) {assert(0 && "not implemented in blas");}
			else                                     {assert(0 && "not implemented in blas");}
		}
		#undef CTXT
	}
	return std::forward<B2D>(b);
//} catch(std::logic_error& le) {
//	using std::to_string;
//	throw std::logic_error{
//		"couldn't do "+std::string(__PRETTY_FUNCTION__)+" of layout a_side="+ static_cast<char>(a_side) +" a_fill="+ static_cast<char>(a_fill) +" a_diag="+ static_cast<char>(a_diag) +" alpha=xx" // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay)
//		+" a_conj="+ to_string(is_conjugated<A2D>{}) +" a_strides="+to_string(stride(a)) +","+ to_string(stride(~a))+" a_sizes="+to_string(size(a)) +","+ to_string(size(~a))
//		+" b_conj="+ to_string(is_conjugated<B2D>{}) +" b_strides="+to_string(stride(b)) +","+ to_string(stride(~b))+" b_sizes="+to_string(size(b)) +","+ to_string(size(~b))
//		+" because " + le.what()
//	};
}

template<class A2D, class B2D>
auto trsm(blas::side a_side, blas::filling a_fill, blas::diagonal a_diag, typename A2D::element_type alpha, A2D const& a, B2D&& b) -> decltype(auto) {
	if constexpr(not is_conjugated<A2D>{}) {return trsm(default_context_of(           a.base() ), a_side, a_fill, a_diag, alpha, a, std::forward<B2D>(b));}
	else                                   {return trsm(default_context_of(underlying(a.base())), a_side, a_fill, a_diag, alpha, a, std::forward<B2D>(b));}
}

template<class Context, class A2D, class B2D>
auto trsm(Context&& ctxt, blas::side a_side, blas::filling a_fill, typename A2D::element_type alpha, A2D const& a, B2D&& b)
->decltype(trsm(std::forward<Context>(ctxt), a_side, a_fill, blas::diagonal::general, alpha, a, std::forward<B2D>(b))) {
	return trsm(std::forward<Context>(ctxt), a_side, a_fill, blas::diagonal::general, alpha, a, std::forward<B2D>(b)); }

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
auto trsm(blas::side a_side, blas::filling a_fill, typename A2D::element_type alpha, A2D const& a, B2D&& b) -> decltype(auto) {
	if constexpr(not is_conjugated<A2D>{}) {return trsm(default_context_of(           a.base() ), a_side, a_fill, alpha, a, std::forward<B2D>(b));}
	else                                   {return trsm(default_context_of(underlying(a.base())), a_side, a_fill, alpha, a, std::forward<B2D>(b));}
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

}  // end namespace boost::multi::blas
