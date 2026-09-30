// Copyright 2020-2025 Alfredo A. Correa

#ifndef BOOST_MULTI_ADAPTORS_LAPACK_GETRF_HPP
#define BOOST_MULTI_ADAPTORS_LAPACK_GETRF_HPP

#include "boost/multi/adaptors/blas/filling.hpp"
#include "boost/multi/adaptors/lapack/core.hpp"

#include "boost/multi/detail/config/NODISCARD.hpp"

#include<cassert>

namespace boost::multi::lapack {

using index = int;

using blas::filling;

template<class Context, class A, class IPIV>
auto getrf(Context&& ctxt, A&& arr, IPIV&& ipiv){
	assert( ipiv.size() == std::min(size(arr), size(~arr)) );
	assert( stride(arr) == 1 );
//  assert( stride(ipiv) == 1 );
	multi::index const i = std::forward<Context>(ctxt).getrf(size(~arr), size(arr), arr.base(), stride(~arr), std::forward<IPIV>(ipiv).data() );
	// if(i == 0) { return arr(); }
	// else       { return arr({0, i - 1}, {0, i - 1}); }
	if(i == 0) { return std::forward<A>(arr)(); }
	return std::forward<A>(arr)({0, i - 1}, {0, i - 1});
}

template<class Context, class LU, class IPIV, class B>
void getrs(Context&& ctxt, LU const& lu, IPIV const& ipiv, B&& barr){
	assert( size(lu) == size(~lu) );
	assert( stride(lu) == 1 );
	assert( size(ipiv) >= size(lu) );
//  assert( stride(ipiv) == 1 );
	assert( stride(barr) == 1 );
	std::forward<Context>(ctxt).getrs('N', size(lu), size(~barr), lu.base(), stride(~lu), ipiv.data(), std::forward<B>(barr).base(), stride(~barr));
}

template<class Context, class LU, class IPIV, class V>
void getrs_one(Context&& ctxt, LU const& lu, IPIV const& ipiv, V&& barr){
	assert( size(lu) == size(~lu) );
	assert( stride(lu) == 1 );
//  assert( stride(ipiv) == 1 );
	assert( stride(barr) == 1 );
	std::forward<Context>(ctxt).getrs('N', size(lu), 1, lu.base(), stride(~lu), ipiv.data(), std::forward<V>(barr).base(), size(lu));
}


template<class A, class IPIV>
auto getrf(A&& arr, IPIV&& ipiv){return getrf(::lapack::context{}, std::forward<A>(arr), std::forward<IPIV>(ipiv));}

template<class LU, class IPIV, class B>
void getrs(LU const& lu, IPIV const& ipiv, B&& barr){return getrs(::lapack::context{}, lu, ipiv, std::forward<B>(barr));}

template<class LU, class IPIV, class B>
void getrs_one(LU const& lu, IPIV const& ipiv, B&& barr){return getrs_one(::lapack::context{}, lu, ipiv, std::forward<B>(barr));}


}  // namespace boost::multi::lapack

#endif


