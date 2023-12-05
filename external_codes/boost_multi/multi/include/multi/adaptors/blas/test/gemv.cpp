// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;-*-
// © Alfredo A. Correa 2020

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi BLAS gemv"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include "config.hpp"

#include "../../../adaptors/blas/gemv.hpp"
#include "../../../array.hpp"

#include "../../../utility.hpp"

#include "../../blas/axpy.hpp"
#include "../../blas/dot.hpp"
#include "../../blas/gemm.hpp"
#include "../../blas/nrm2.hpp"

#include<random>

namespace multi = boost::multi;
namespace blas = multi::blas;

template<class T> void what(T&&) = delete;

template<class M, class VI, class VO>
void MV(M const& a, VI const& x, VO&& y) {  // NOLINT(readability-identifier-naming,readability-identifier-length) BLAS naming
	std::transform(
		begin(a), end(a), begin(y),
		[&x](auto&& row){return std::inner_product(begin(row), end(row), begin(x), 0.);}
	);
}

BOOST_AUTO_TEST_CASE(multi_blas_gemv) {
	multi::array<double, 2> const a = {  // NOLINT(readability-identifier-length) BLAS naming
		{ 9., 24., 30., 9.},
		{ 4., 10., 12., 7.},
		{14., 16., 36., 1.}
	};
	multi::array<double, 1> const x = {1.1, 2.1, 3.1, 4.1};  // NOLINT(readability-identifier-length) BLAS naming
	{
		multi::array<double, 1>       y(multi::extensions_t<1>{multi::iextension{size(a)}});  // NOLINT(readability-identifier-length) BLAS naming
		blas::gemv_n(1., begin(a), size(a), begin(x), 0., begin(y));
		BOOST_REQUIRE_CLOSE( y[1] , 91.3                , 0.0001 );
		BOOST_REQUIRE_CLOSE( y[2] , +blas::dot(a[2], x) , 0.0001 );
	}
	{
		multi::array<double, 1>       y(multi::extensions_t<1>{multi::iextension{size(a)}});  // NOLINT(readability-identifier-length) BLAS naming
		multi::array<double, 2> const aT = ~a;
		blas::gemv_n(1., begin(~aT), size(~aT), begin(x), 0., begin(y));
		BOOST_REQUIRE_CLOSE( y[1] , 91.3               , 0.0001 );
		BOOST_REQUIRE_CLOSE( y[2] , +blas::dot(a[2], x), 0.0001 );
	}
	{
		multi::array<double, 1>       y(multi::extensions_t<1>{multi::iextension{size(a)}});  // NOLINT(readability-identifier-length) BLAS naming
		auto mv = blas::gemv(1., a, x);
		copy_n(mv.begin(), mv.size(), y.begin());
		BOOST_REQUIRE_CLOSE( y[1] , 91.3 , 0.00001 );

		multi::array<double, 1> w2(multi::extensions_t<1>{multi::iextension{size(a)}});
		MV(a, x, w2);
		BOOST_REQUIRE_CLOSE( w2[0] , y[0], 0.00001 );
	}
	{
		multi::array<double, 1>       y(multi::extensions_t<1>{multi::iextension{size(a)}});  // NOLINT(readability-identifier-length) BLAS naming
		y = blas::gemv(1., a, x);
		BOOST_REQUIRE_CLOSE( y[1] , 91.3 , 0.00001 );
	}
	{
		multi::array<double, 1> y = blas::gemv(1., a, x);  // NOLINT(readability-identifier-length) BLAS naming
		BOOST_REQUIRE_CLOSE( y[1] , 91.3 , 0.00001 );
	}
	{
		multi::array<double, 1>       y(multi::extensions_t<1>{multi::iextension{size(a)}}, 0.);  // NOLINT(readability-identifier-length) BLAS naming
		y += blas::gemv(1., a, x);
		BOOST_REQUIRE_CLOSE( y[1] , 91.3 , 0.00001 );
	}
	 {
		multi::array<double, 1> y = {4., 5., 6.};  // NOLINT(readability-identifier-length) BLAS naming
		blas::gemv(1.1, a, x, 1., y);              // y = a*M*x + b*y
		BOOST_REQUIRE_CLOSE( y[1] , 105.43 , 0.00001 );
	}
}

BOOST_AUTO_TEST_CASE(multi_blas_gemv_real) {
	namespace blas = multi::blas;

	using std::abs;
	multi::array<double, 2> const a = {  // NOLINT(readability-identifier-length) BLAS naming
		{ 9., 24., 30., 9.},
		{ 4., 10., 12., 7.},
		{14., 16., 36., 1.}
	};
	multi::array<double, 1> const x = {1.1, 2.1, 3.1, 4.1};  // NOLINT(readability-identifier-length) BLAS naming
	{
		multi::array<double, 1> y = {4., 5., 6.};  // NOLINT(readability-identifier-length) BLAS naming
		double const alpha = 1.1;
		double const beta = 1.2;
		blas::gemv(alpha, a, x, beta, y);  // y = a*M*x + b*y

		multi::array<double, 1> const y3 = {214.02, 106.43, 188.37};
		BOOST_REQUIRE( abs(y[1] - y3[1]) < 2e-14 );
	}
	{
		auto Y = +blas::gemv(1., a, x);  // NOLINT(readability-identifier-length) BLAS naming
		BOOST_REQUIRE_CLOSE( Y[0] , +blas::dot(a[0], x) , 0.00001 );
		BOOST_REQUIRE_CLOSE( Y[1] , +blas::dot(a[1], x) , 0.00001 );
		BOOST_REQUIRE_CLOSE( Y[2] , +blas::dot(a[2], x) , 0.00001 );
	}
	{
		multi::array<double, 1> const x = {1., 2., 3.};  // NOLINT(readability-identifier-length) BLAS naming
		multi::array<double, 1> const y = {4., 5., 6.};  // NOLINT(readability-identifier-length) BLAS naming
		multi::array<double, 1> const dot = blas::gemv(1., multi::array<double, 2>({x}), y);
		BOOST_REQUIRE( dot[0] == blas::dot(x, y) );
	}
	{
		using blas::operators::operator%;
		using blas::operators::operator-;
		using blas::operators::operator^;
		BOOST_REQUIRE_SMALL( ((~+~a)%x - a%x)^2 , 1e-13 );
	}
}

BOOST_AUTO_TEST_CASE(multi_blas_gemv_real_complex) {
	namespace blas = multi::blas;
	using complex = std::complex<double>; //#define I *std::complex<double>(0, 1)
	using std::abs;
	multi::array<complex, 2> const M = {  // NOLINT(readability-identifier-length) BLAS naming
		{ 9., 24., 30., 9.},
		{ 4., 10., 12., 7.},
		{14., 16., 36., 1.}
	};
	multi::array<complex, 1> const X = {1.1, 2.1, 3.1, 4.1};  // NOLINT(readability-identifier-length) BLAS naming
	{
		multi::array<complex, 1> Y = {4., 5., 6.};  // NOLINT(readability-identifier-length) BLAS naming
		double const alpha = 1.1;
		double const beta = 1.2;
		blas::gemv(alpha, M, X, beta, Y); // y = a*M*x + b*y

		multi::array<complex, 1> const Y3 = {214.02, 106.43, 188.37};

		using blas::operators::operator-;
		double const n2{blas::nrm2(Y - Y3)};
		BOOST_REQUIRE_SMALL( n2 , 1e-13);
	}
}

#if CUDA_FOUND
#include<thrust/complex.h>
BOOST_AUTO_TEST_CASE(multi_blas_gemv_real_complex_thrust) {
	namespace blas = multi::blas;
	using complex = thrust::complex<double>; //#define I *std::complex<double>(0, 1)
	using std::abs;
	multi::array<complex, 2> const M = {
		{ 9., 24., 30., 9.},
		{ 4., 10., 12., 7.},
		{14., 16., 36., 1.}
	};
	multi::array<complex, 1> const X = {1.1, 2.1, 3.1, 4.1};
	{
		multi::array<complex, 1> Y = {4., 5., 6.};
		double const a = 1.1;
		double const b = 1.2;
		blas::gemv(a, M, X, b, Y); // y = a*M*x + b*y

		multi::array<complex, 1> const Y3 = {214.02, 106.43, 188.37};
	}
	 {
		multi::array<complex, 1> Y = {4., 5., 6.};
		blas::gemv(1.1, M, X, 1., Y); // y = a*M*x + b*y
		BOOST_REQUIRE( Y[1] == 105.43 );
	}
}
#endif

BOOST_AUTO_TEST_CASE(multi_blas_gemv_complex) {
	namespace blas = multi::blas;
	using complex = std::complex<double>; std::complex<double> const I{0, 1};  // NOLINT(readability-identifier-length) imag unit

	using std::abs;
	multi::array<complex, 2> const a = {  // NOLINT(readability-identifier-length) BLAS naming
		{2. + 3.*I, 2. + 1.*I, 1. + 2.*I},
		{4. + 2.*I, 2. + 4.*I, 3. + 1.*I},
		{7. + 1.*I, 1. + 5.*I, 0. + 3.*I}
	};
	multi::array<complex, 1> const x = {1. + 2.*I, 2. + 1.*I, 9. + 2.*I};  // NOLINT(readability-identifier-length) BLAS naming
	BOOST_REQUIRE(( +blas::gemv(1., a, x) == multi::array<complex, 1>{4. + 31.*I, 25. + 35.*I, -4. + 53.*I} ));

	auto aT = +~a;
	BOOST_REQUIRE(( +blas::gemv(1., ~aT, x) == multi::array<complex, 1>{4. + 31.*I, 25. + 35.*I, -4. + 53.*I} ));

	BOOST_REQUIRE( +blas::gemv(1., ~a, x) == (multi::array<complex, 1>{63. + 38.*I, -1. + 62.*I, -4. + 36.*I}) );
	BOOST_REQUIRE( +blas::gemv(1., ~a, x) == +blas::gemv(1., aT, x) );
}

BOOST_AUTO_TEST_CASE(multi_blas_gemv_temporary) {
	using complex = std::complex<double>;

	multi::array<complex, 2> const A = {  // NOLINT(readability-identifier-length) BLAS naming
		{1., 0., 0.},
		{0., 1., 0.},
		{0., 0., 1.}
	};

	auto const B = [](auto array) {  // NOLINT(readability-identifier-length) BLAS naming
		auto rand = [gauss = std::normal_distribution<>{}, gen = std::mt19937{1}]() mutable {return complex{gauss(gen), gauss(gen)};};  // NOLINT(cert-msc32-c,cert-msc51-cpp) test purposes
		std::generate(array.elements().begin(), array.elements().end(), rand);
		return array;
	}(multi::array<complex, 2>({3, 3}));

	using blas::operators::operator*;
	using blas::operators::operator-;
	using blas::operators::operator^;
	BOOST_REQUIRE( (((A*B)[0] - B[0])^2) == 0. );
	BOOST_REQUIRE( (((A*B)[1] - B[1])^2) == 0. );
	BOOST_REQUIRE( (((A*B)[2] - B[2])^2) == 0. );
}

BOOST_AUTO_TEST_CASE(multi_blas_gemv_context) {
	multi::array<double, 2> const a = {  // NOLINT(readability-identifier-length) BLAS naming
		{ 9., 24., 30., 9.},
		{ 4., 10., 12., 7.},
		{14., 16., 36., 1.}
	};
	multi::array<double, 1> const x = {1.1, 2.1, 3.1, 4.1};  // NOLINT(readability-identifier-length) conventional name in BLAS

	blas::context ctxt;
	{
		multi::array<double, 1>       y(multi::extensions_t<1>{multi::iextension{size(a)}});  // NOLINT(readability-identifier-length) conventional name in BLAS
		blas::gemv_n(ctxt, 1., begin(a), size(a), begin(x), 0., begin(y));
		BOOST_REQUIRE_CLOSE( y[1] , 91.3 , 0.0001 );
		BOOST_REQUIRE_CLOSE( y[2] , +blas::dot(a[2], x) , 0.0001 );
	}
	{
		multi::array<double, 1>       y(multi::extensions_t<1>{multi::iextension{size(a)}});  // NOLINT(readability-identifier-length) conventional name in BLAS
		multi::array<double, 2> const aT = ~a;
		blas::gemv_n(ctxt, 1., begin(~aT), size(~aT), begin(x), 0., begin(y));
		BOOST_REQUIRE_CLOSE( y[1] , 91.3 , 0.00001 );
		BOOST_REQUIRE_CLOSE( y[2] , +blas::dot(a[2], x) , 0.00001 );
	}
	{
		multi::array<double, 1>       y(multi::extensions_t<1>{multi::iextension{size(a)}});  // NOLINT(readability-identifier-length) conventional name in BLAS
		auto&& mv = blas::gemv(ctxt, 1., a, x);
		copy_n(mv.begin(), mv.size(), y.begin());
		BOOST_REQUIRE_CLOSE( y[1] , 91.3 , 0.00001 );
	}
	{
		multi::array<double, 1>       y(multi::extensions_t<1>{multi::iextension{size(a)}});  // NOLINT(readability-identifier-length) conventional name in BLAS
		y = blas::gemv(ctxt, 1., a, x);
		BOOST_REQUIRE_CLOSE( y[1] , 91.3 , 0.00001 );
	}
	{
		multi::array<double, 1> y = blas::gemv(ctxt, 1., a, x);  // NOLINT(readability-identifier-length) conventional name in BLAS
		BOOST_REQUIRE_CLOSE( y[1] , 91.3 , 0.00001 );
	}
	{
		multi::array<double, 1>       y(multi::extensions_t<1>{multi::iextension{size(a)}}, 0.);  // NOLINT(readability-identifier-length) conventional name in BLAS
		y += blas::gemv(ctxt, 1., a, x);
		BOOST_REQUIRE_CLOSE( y[1] , 91.3, 0.00001 );
	}
	 {
		multi::array<double, 1> y = {4., 5., 6.};  // NOLINT(readability-identifier-length) conventional name in BLAS
		y += blas::gemv(ctxt, 1.1, a, x);
		BOOST_REQUIRE_CLOSE( y[1] , 105.43, 0.00001 );
	}
}
