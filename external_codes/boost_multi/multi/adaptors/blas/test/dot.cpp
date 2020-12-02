#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;-*-
$CXX $0 -o $0x `pkg-config --libs blas` -lcudart -lcublas -lboost_unit_test_framework&&$0x&&rm $0x;exit
#endif
// © Alfredo A. Correa 2019-2020

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi cuBLAS dot"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include "../../blas/dot.hpp"

#include "../../../array.hpp"
#include "../../blas/cuda.hpp"

#include "../../../adaptors/cuda.hpp"
#include "../../../complex.hpp"

#include<type_traits>
#include<thrust/complex.h>
namespace std{
	template<> struct is_trivially_constructible<thrust::complex<double>, thrust::complex<double>> : std::true_type{};
	template<> struct is_trivially_assignable<thrust::complex<double>&, thrust::complex<double>> : std::true_type{};
}

#include<complex>
#include<cassert>
#include<numeric>

using std::cout;
namespace multi = boost::multi;
namespace blas = multi::blas;

BOOST_AUTO_TEST_CASE(blas_dot){
	namespace cuda = multi::cuda;
	{
		multi::array<double, 2> CA = {
			{1.,  2.,  3.,  4.},
			{5.,  6.,  7.,  8.},
			{9., 10., 11., 12.}
		};
		auto d = blas::dot(CA[1], CA[2]);
		BOOST_REQUIRE( d() == std::inner_product(begin(CA[1]), begin(CA[2]), end(CA[1]), 0.) );
	}
	{
		multi::array<float, 1> const A = {1.,2.,3.};
		multi::array<float, 1> const B = {1.,2.,3.};
		auto f = blas::dot(A, B); // sdot
		BOOST_REQUIRE( f() == std::inner_product(begin(A), end(A), begin(B), float{0}) );
	}
	using complex = std::complex<double>; complex const I{0, 1};
	{
		multi::array<complex, 1> const A = {I, 2.*I, 3.*I};
		BOOST_REQUIRE( blas::dot(A, A)() == std::inner_product(begin(A), end(A), begin(A), std::complex<double>(0)) );
	}
	{
		multi::array<complex, 1> const A = {I, 2.*I, 3.*I};
		BOOST_REQUIRE( blas::dot(A, blas::C(A))() == std::inner_product(begin(A), end(A), begin(A), std::complex<double>(0), std::plus<>{}, [](auto&& a, auto&& b){return a*conj(b);}) );
	}
	{
		multi::array<complex, 1> const a = {1. + I, 2. + 3.*I, 3. + 2.*I, 4. - 9.*I};
		multi::array<complex, 1> const b = {5. + 2.*I, 6. + 6.*I, 7. + 2.*I, 8. - 3.*I};
		{
			multi::array<complex, 0> c;
			blas::dot(a, b, c);
			BOOST_REQUIRE( c() == 19. - 27.*I );
		}
	}
	{
		cuda::array<complex, 1> const acu = {1. + I, 2. + 3.*I, 3. + 2.*I, 4. - 9.*I};
		cuda::array<complex, 1> const bcu = {5. + 2.*I, 6. + 6.*I, 7. + 2.*I, 8. - 3.*I};

		{
			cuda::array<complex, 0> ccu;
			blas::dot(acu, bcu, ccu);
			BOOST_REQUIRE( ccu() == 19. - 27.*I );
		}
		BOOST_REQUIRE( blas::C(bcu)[1] == 6. - 6.*I );
		{
			cuda::array<complex, 0> ccu;
			static_assert( multi::blas::is_complex_array<multi::array<complex, 1>>{}, "!" );
			static_assert( multi::blas::is_complex_array<cuda::array<complex, 1>>{}, "!" );
			blas::dot(acu, blas::C(bcu), ccu);
			BOOST_REQUIRE( ccu() == 121. - 43.*I );
		}
		{
			auto const ccu = blas::dot(acu, blas::C(bcu));
			BOOST_REQUIRE( ccu() == 121. - 43.*I );
		}
		{
			cuda::array<complex, 1> ccu = {1, 2, 3};
			blas::dot(acu, blas::C(bcu), ccu[0]);
			BOOST_REQUIRE( ccu[0] == 121. - 43.*I );
		}
		{
			cuda::array<complex, 2> ccu({1, 1});
			blas::dot(acu, blas::C(bcu), ccu[0][0]);
			BOOST_REQUIRE( ccu[0][0] == 121. - 43.*I );
		}
	}
	{
		namespace cuda = multi::cuda;
		cuda::managed::array<complex, 1> const amcu = {1. + I, 2. + 3.*I, 3. + 2.*I, 4. - 9.*I};
		cuda::managed::array<complex, 1> const bmcu = {5. + 2.*I, 6. + 6.*I, 7. + 2.*I, 8. - 3.*I};
		{
			cuda::managed::array<complex, 0> cmcu;
			blas::dot(amcu, bmcu, cmcu);
			BOOST_REQUIRE( cmcu() == 19.- I*27. );
		}
		{
			cuda::array<complex, 1> cmcu = {1, 2, 3};
			blas::dot(amcu, blas::C(bmcu), cmcu[0]);
			BOOST_REQUIRE( cmcu[0] == complex(121., -43.) );
		}
	}
	{
		using complex = std::complex<double>; complex const I{0, 1};
		cuda::array<complex, 1> const acu = {1. + I, 2. + 3.*I, 3. + 2.*I, 4. - 9.*I};
		cuda::array<complex, 1> const bcu = {5. + 2.*I, 6. + 6.*I, 7. + 2.*I, 8. - 3.*I};
		{
			cuda::array<complex, 0> ccu;
			blas::dot(acu, bcu, ccu);
			BOOST_REQUIRE( ccu() == 19. - 27.*I );
		}
	}
	{
		using complex = thrust::complex<double>; complex const I{0, 1};
		cuda::managed::array<complex, 1> const acu = {1. + I, 2. + 3.*I, 3. + 2.*I, 4. - 9.*I};
		cuda::managed::array<complex, 1> const bcu = {5. + 2.*I, 6. + 6.*I, 7. + 2.*I, 8. - 3.*I};
		{
			cuda::managed::array<complex, 0> ccu;
			blas::dot(acu, bcu, ccu);
			BOOST_REQUIRE( ccu() == 19. - 27.*I );
		}
	}
}

#if 0
BOOST_AUTO_TEST_CASE(multi_blas_dot_impl_thrust_complex){
	namespace blas = multi::blas;

	using complex = thrust::complex<double>; complex const I{0, 1};
	multi::array<complex, 2> const A = {
		{1. +    I,  2. + 3.*I,  3.+2.*I,  4.-9.*I},
		{5. + 2.*I,  6. + 6.*I,  7.+2.*I,  8.-3.*I},
		{9. + 1.*I, 10. + 9.*I, 11.+1.*I, 12.+2.*I}
	};
	{
	//	complex c; blas::dot(A[1], A[2], c);
	//	BOOST_REQUIRE( c == std::inner_product(begin(A[1]), end(A[1]), begin(A[2]), complex{0}) );
	}
}
#endif

