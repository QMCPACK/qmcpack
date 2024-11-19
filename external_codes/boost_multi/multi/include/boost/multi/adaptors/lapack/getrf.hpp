// Copyright 2020-2024 Alfredo A. Correa

#ifndef BOOST_MULTI_ADAPTORS_LAPACK_GEQRF_HPP
#define BOOST_MULTI_ADAPTORS_LAPACK_GEQRF_HPP

#include <boost/multi/adaptors/lapack/core.hpp>
#include <boost/multi/adaptors/blas/filling.hpp>

#include <boost/multi/detail/config/NODISCARD.hpp"

#include<cassert>

//#include<lapacke.h>

namespace boost{namespace multi{namespace lapack{

using index = int;

using blas::filling;

template<class Context, class A, class IPIV>
auto getrf(Context&& ctxt, A&& a, IPIV&& ipiv){
	assert( ipiv.size() == std::min(size(a), size(~a)) );
	assert( stride(a) == 1 );
//  assert( stride(ipiv) == 1 );
	multi::index i = std::forward<Context>(ctxt).getrf(size(~a), size(a), base(a), stride(~a), ipiv.data() );
	if(i == 0) return a();
	else       return a({0, i - 1}, {0, i - 1});
}

template<class Context, class LU, class IPIV, class B>
void getrs(Context&& ctxt, LU const& lu, IPIV const& ipiv, B&& b){
	assert( size(lu) == size(~lu) );
	assert( stride(lu) == 1 );
	assert( size(ipiv) >= size(lu) );
//  assert( stride(ipiv) == 1 );
	assert( stride(b) == 1 );
	std::forward<Context>(ctxt).getrs('N', size(lu), size(~b), base(lu), stride(~lu), ipiv.data(), base(b), stride(~b));
}

template<class Context, class LU, class IPIV, class V>
void getrs_one(Context&& ctxt, LU const& lu, IPIV const& ipiv, V&& b){
	assert( size(lu) == size(~lu) );
	assert( stride(lu) == 1 );
//  assert( stride(ipiv) == 1 );
	assert( stride(b) == 1 );
	std::forward<Context>(ctxt).getrs('N', size(lu), 1, base(lu), stride(~lu), ipiv.data(), base(b), size(lu));
}


template<class A, class IPIV>
auto getrf(A&& a, IPIV&& ipiv){return getrf(::lapack::context{}, std::forward<A>(a), std::forward<IPIV>(ipiv));}

template<class LU, class IPIV, class B>
void getrs(LU const& lu, IPIV const& ipiv, B&& b){return getrs(::lapack::context{}, lu, ipiv, std::forward<B>(b));}

template<class LU, class IPIV, class B>
void getrs_one(LU const& lu, IPIV const& ipiv, B&& b){return getrs_one(::lapack::context{}, lu, ipiv, std::forward<B>(b));}


}}}

#endif


