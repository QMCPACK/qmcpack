#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;-*-
$CXX $0 -o $0x `pkg-config --libs blas` -lcudart -lcublas -lboost_unit_test_framework&&$0x&&rm $0x;exit
#endif
// © Alfredo A. Correa 2019-2020

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi BLAS operations and cuda"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include "../../blas/dot.hpp"

#include "../../../array.hpp"
#include "../../blas/cuda.hpp"

#include "../../../adaptors/cuda.hpp"
#include "../../../complex.hpp"

#include<complex>
#include<cassert>
#include<numeric>

using std::cout;
namespace multi = boost::multi;
namespace blas = multi::blas;

using complex = std::complex<double>; constexpr complex I{0, 1};

BOOST_AUTO_TEST_CASE(blas_conjugated_cpu){
	multi::array<complex, 1> const a = {5. + 2.*I, 6. + 6.*I, 7. + 2.*I, 8. - 3.*I};
	BOOST_REQUIRE( blas::C(a)[1] == conj(a[1]) );

	namespace cuda = multi::cuda;

	cuda::array<complex, 1> const agpu = {5. + 2.*I, 6. + 6.*I, 7. + 2.*I, 8. - 3.*I};
	BOOST_REQUIRE( blas::C(agpu)[1] == conj(agpu[1]) );
}

BOOST_AUTO_TEST_CASE(blas_conjugated_gpu){
#if 0
	cuda::array<complex, 1> const acu = {1. +    I, 2. + 3.*I, 3. + 2.*I, 4. - 9.*I};
	cuda::array<complex, 1> const bcu = {5. + 2.*I, 6. + 6.*I, 7. + 2.*I, 8. - 3.*I};

	{
		cuda::array<complex, 0> ccu;
		blas::dot(acu, bcu, ccu);
		BOOST_REQUIRE( ccu() == 19. - 27.*I );
	}
	BOOST_REQUIRE( blas::C(bcu)[1] == 2. - 3.*I );
#endif
}

