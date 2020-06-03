#ifdef COMPILATION_INSTRUCTIONS
$CXX -Wall -Wextra -Wpedantic `#-Wfatal-errors` $0 -o $0x `pkg-config --libs blas` -lcudart -lcublas -lboost_unit_test_framework&&$0x&&rm $0x;exit
#endif
// © Alfredo A. Correa 2019-2020

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi cuBLAS dot"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include "../../blas/dot.hpp"

#include "../../../array.hpp"
#include "../../blas/cuda.hpp"

#include "../../../adaptors/cuda.hpp"

#include<complex>
#include<cassert>
#include<numeric>

using std::cout;
namespace multi = boost::multi;
namespace blas = multi::blas;

BOOST_AUTO_TEST_CASE(blas_dot){
	{
		multi::array<double, 2> CA = {
			{1.,  2.,  3.,  4.},
			{5.,  6.,  7.,  8.},
			{9., 10., 11., 12.}
		};
		using blas::dot;
		auto d = dot(CA[1], CA[2]);
		assert( d() == std::inner_product(begin(CA[1]), begin(CA[2]), end(CA[1]), 0.) );
	}
	{
		multi::array<float, 1> const A = {1.,2.,3.};
		multi::array<float, 1> const B = {1.,2.,3.};
		using blas::dot;
		auto f = dot(A, B); // sdot
		assert( f() == std::inner_product(begin(A), end(A), begin(B), float{0}) );
	}
	using complex = std::complex<double>; complex const I{0, 1};
	{
		multi::array<complex, 1> const A = {I, 2.*I, 3.*I};
		using blas::dot;
		assert( dot(A, A)() == std::inner_product(begin(A), end(A), begin(A), std::complex<double>(0)) );
	}
	{
		multi::array<complex, 1> const A = {I, 2.*I, 3.*I};
		using blas::dot;
		using blas::conjugated;
		assert( dot(A, conjugated(A))() == std::inner_product(begin(A), end(A), begin(A), std::complex<double>(0), std::plus<>{}, [](auto&& a, auto&& b){return a*conj(b);}) );
	}
	{
		multi::array<complex, 1> const a = {1. + I, 2. + 3.*I, 3. + 2.*I, 4. - 9.*I};
		multi::array<complex, 1> const b = {5. + 2.*I, 6. + 6.*I, 7. + 2.*I, 8. - 3.*I};

		using blas::conjugated;
		using blas::dot;
		{
			multi::array<complex, 0> c;
			dot(a, b, c);
			BOOST_REQUIRE( c() == 19. - 27.*I );
		}
	}
	{
		namespace cuda = multi::cuda;
		cuda::array<complex, 1> const acu = {1. + I, 2. + 3.*I, 3. + 2.*I, 4. - 9.*I};
		cuda::array<complex, 1> const bcu = {5. + 2.*I, 6. + 6.*I, 7. + 2.*I, 8. - 3.*I};

		using blas::conjugated;
		{
			cuda::array<complex, 0> ccu;
			blas::dot(acu, bcu, ccu);
			BOOST_REQUIRE( ccu() == 19. - 27.*I );
		}
		{
			cuda::array<complex, 0> ccu;
			blas::dot(acu, conjugated(bcu), ccu);
			BOOST_REQUIRE( ccu() == 121. - 43.*I );
		}
		{
			auto const ccu = blas::dot(acu, conjugated(bcu));
			BOOST_REQUIRE( ccu() == 121. - 43.*I );
		}
		{
			cuda::array<complex, 1> ccu = {1, 2, 3};
			blas::dot(acu, conjugated(bcu), ccu[0]);
			BOOST_REQUIRE( ccu[0] == 121. - 43.*I );
		}
		{
			cuda::array<complex, 2> ccu({1, 1});
			blas::dot(acu, conjugated(bcu), ccu[0][0]);
			BOOST_REQUIRE( ccu[0][0] == 121. - 43.*I );
		}
	}
	{
		namespace cuda = multi::cuda;
		cuda::managed::array<complex, 1> const amcu = {1. + I, 2. + 3.*I, 3. + 2.*I, 4. - 9.*I};
		cuda::managed::array<complex, 1> const bmcu = {5. + 2.*I, 6. + 6.*I, 7. + 2.*I, 8. - 3.*I};

		using blas::conjugated;
		{
			cuda::managed::array<complex, 0> cmcu;
			blas::dot(amcu, bmcu, cmcu);
			BOOST_REQUIRE( cmcu() == 19.- I*27. );
		}
	}
#if 0
	{
		{
			cuda::array<complex, 1> ccu = {1, 2, 3};
			dot(acu, conjugated(bcu), ccu[0]);
			BOOST_REQUIRE( ccu[0] == complex(227, 0) );
		}
	}
#endif
}

