// Copyright 2023-2025 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/multi/adaptors/blas/herk.hpp>
#include <boost/multi/adaptors/thrust.hpp>
#include <thrust/complex.h>

#include <boost/core/lightweight_test.hpp>
#define BOOST_AUTO_TEST_CASE(CasenamE) /**/

namespace multi = boost::multi;
namespace blas  = multi::blas;

using complex = thrust::complex<double>;
complex const I{0, 1};

int main() {
	BOOST_AUTO_TEST_CASE(multi_blas_herk){
		multi::array<thrust::complex<double>, 2> const a = {
			{1.0 + 3.0 * I, 3.0 - 2.0 * I, 4.0 + 1.0 * I},
			{9.0 + 1.0 * I, 7.0 - 8.0 * I, 1.0 - 3.0 * I}
		};
		multi::thrust::cuda::array<thrust::complex<double>, 2> const a_gpu = a;

		{
			multi::array<thrust::complex<double>, 2> c({2, 2}, 9999.0);
			blas::herk(1.0, a, c);
			BOOST_TEST( c[1][0] == complex(50.0, -49.0) );
			BOOST_TEST( c[0][1] == complex(50.0, +49.0) );

			multi::array<thrust::complex<double>, 2> const c_copy = blas::herk(1.0, a);
			BOOST_TEST( c == c_copy );
		}
		{
			multi::array<thrust::complex<double>, 2> c({3, 3}, 9999.0);
			blas::herk(1.0, blas::H(a), c);
			BOOST_TEST( c[2][1] == complex(41, +2) );
			BOOST_TEST( c[1][2] == complex(41, -2) );

			multi::array<thrust::complex<double>, 2> const c_copy = blas::herk(1., blas::H(a));
			BOOST_TEST( c_copy == c );
		}
	}

	return boost::report_errors();
}
