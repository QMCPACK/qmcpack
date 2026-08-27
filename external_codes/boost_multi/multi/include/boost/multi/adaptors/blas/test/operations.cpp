// Copyright 2019-2026 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/core/lightweight_test.hpp>

#include "../../blas/dot.hpp"

#include "../../../array.hpp"
#include "../../blas/cuda.hpp"

#include "../../../adaptors/cuda.hpp"
#include "../../../complex.hpp"

#include<complex>
#include<cassert>
#include<numeric>

namespace multi = boost::multi;
namespace blas = multi::blas;

using complex = std::complex<double>; constexpr complex I{0.0, 1.0};

auto main() -> int {  // NOLINT(bugprone-exception-escape)
	// BOOST_AUTO_TEST_CASE(blas_conjugated_cpu)
	{
		multi::array<complex, 1> const a = {5.0 + 2.0*I, 6.0 + 6.0*I, 7.0 + 2.0*I, 8.0 - 3.0*I};
		BOOST_TEST( blas::C(a)[1] == conj(a[1]) );

		namespace cuda = multi::cuda;

		cuda::array<complex, 1> const agpu = {5.0 + 2.0*I, 6.0 + 6.0*I, 7.0 + 2.0*I, 8.0 - 3.0*I};
		BOOST_TEST( blas::C(agpu)[1] == conj(agpu[1]) );
	}

	// BOOST_AUTO_TEST_CASE(blas_conjugated_gpu)
	// {
	// 	cuda::array<complex, 1> const acu = {1.0 +     I, 2.0 + 3.0*I, 3.0 + 2.0*I, 4.0 - 9.0*I};
	// 	cuda::array<complex, 1> const bcu = {5.0 + 2.0*I, 6.0 + 6.0*I, 7.0 + 2.0*I, 8.0 - 3.0*I};

	// 	{
	// 		cuda::array<complex, 0> ccu;
	// 		blas::dot(acu, bcu, ccu);
	// 		BOOST_TEST( ccu() == 19.0 - 27.0*I );
	// 	}
	// 	BOOST_TEST( blas::C(bcu)[1] == 2.0 - 3.0*I );
	// }

	return boost::report_errors();
}
