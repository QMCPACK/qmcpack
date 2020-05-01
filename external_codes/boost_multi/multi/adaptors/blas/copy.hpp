#ifdef COMPILATION_INSTRUCTIONS
(echo '#include"'$0'"'>$0.cpp)&&$CXX -D_TEST_MULTI_ADAPTORS_BLAS_COPY $0.cpp -o $0x `pkg-config --cflags --libs blas` -lboost_unit_test_framework&&$0x&&rm $0x $0.cpp;exit
#endif
// © Alfredo A. Correa 2019-2020
#ifndef MULTI_ADAPTORS_BLAS_COPY_HPP
#define MULTI_ADAPTORS_BLAS_COPY_HPP

#include "../blas/core.hpp"
#include "../../config/NODISCARD.hpp"

namespace boost{
namespace multi{
namespace blas{

using core::copy;

template<class It, typename Size, class OutIt>
auto copy_n(It first, Size n, OutIt d_first)
->decltype(copy(n, base(first), stride(first), base(d_first), stride(d_first)), d_first + n){
	return copy(n, base(first), stride(first), base(d_first), stride(d_first)), d_first + n;}

template<class X1D, class Y1D>
Y1D&& copy(X1D const& x, Y1D&& y){assert(size(x)==size(y)); assert(offset(x)==0 and offset(y)==0);
	copy(size(x), base(x), stride(x), base(y), stride(y));
	return std::forward<Y1D>(y);
}

template<class X1D, class Ret = typename X1D::decay_type> // TODO multi::array_traits<X1D>::decay_type
NODISCARD("a copied matrix should be assigned")
Ret copy(X1D const& x){
	assert( not offset(x) );
	return copy(x, Ret(size(x), get_allocator(x)));
}

}}}

#if _TEST_MULTI_ADAPTORS_BLAS_COPY

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi cuBLAS gemm"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include "../../array.hpp"
#include "../../utility.hpp"

#include<complex>

namespace multi = boost::multi;
namespace blas = multi::blas;

using complex = std::complex<double>;
complex const I{0, 1};

BOOST_AUTO_TEST_CASE(multi_blas_copy){
	{
		multi::array<double, 1> const A = {1., 2., 3., 4.};
		multi::array<double, 1> B = {5., 6., 7., 8.};
		blas::copy(A, B);
		BOOST_REQUIRE( B == A );
	}
	{
		using complex = std::complex<double>;
		multi::array<complex, 1> const A = {1., 2., 3., 4.};
		multi::array<complex, 1> B = {5., 6., 7., 8.};
		blas::copy(A, B);
		BOOST_REQUIRE( B == A );		
	}
	{
		multi::array<double, 2> const A = {
			{1., 2., 3.},
			{4., 5., 6.},
			{7., 8., 9.}
		};
		multi::array<double, 1> B(3);
		blas::copy(rotated(A)[0], B);
		BOOST_REQUIRE( B == rotated(A)[0] );
	}
	{
		multi::array<complex, 2> const A = {
			{1., 2., 3.},
			{4., 5., 6.},
			{7., 8., 9.}
		};
		multi::array<complex, 1> B = blas::copy(rotated(A)[2]);
		BOOST_REQUIRE( size(B) == 3 );
		BOOST_REQUIRE( B[0] == 3. );
	}
}

#endif
#endif

