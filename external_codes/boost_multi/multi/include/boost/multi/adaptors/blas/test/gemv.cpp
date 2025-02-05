// Copyright 2020-2024 Alfredo A. Correa

#include <boost/test/unit_test.hpp>

#include <boost/mpl/list.hpp>

#include "../../../adaptors/blas/gemv.hpp"
#include "../../../array.hpp"

#include "../../../utility.hpp"

#include "../../blas/axpy.hpp"
#include "../../blas/dot.hpp"
#include "../../blas/gemm.hpp"
#include "../../blas/nrm2.hpp"

#include <random>

namespace multi = boost::multi;
namespace blas  = multi::blas;

using fp_types = boost::mpl::list<double, float>;  // old versions of Boost.Test need MPL Type lists explicitly

template<class M, class VI, class VO>
auto MV(M const& a, VI const& x, VO&& y) -> VO&& {  // NOLINT(readability-identifier-naming,readability-identifier-length) BLAS naming
	std::transform(
		begin(a), end(a), begin(y),
		[&x](auto const& row) { return std::inner_product(begin(row), end(row), begin(x), 0.0); }
	);
	return std::forward<VO>(y);
}

// #ifdef _MULTI_USING_BLAS_MKL
// #include <mkl/mkl_service.h>  // for mkl_free_buffers
// struct Fixture {
//   Fixture()   {mkl_disable_fast_mm(); }  // this is reported to solve memory leaks, but it doesn't with BLA_VENDOR=Intel10_64ilp (non seq) and INTEL_MKL_VERSION 20200004
//   ~Fixture()  { mkl_free_buffers(); }  // this is reported to solve memory leaks, but it doesn't with BLA_VENDOR=Intel10_64ilp (non seq) and INTEL_MKL_VERSION 20200004
// };

// BOOST_GLOBAL_FIXTURE(Fixture);
// #endif

BOOST_AUTO_TEST_CASE_TEMPLATE(multi_blas_gemv, T, fp_types) {
	// NOLINTNEXTLINE(readability-identifier-length) BLAS naming
	multi::array<T, 2> const a = {
		{ 9.0, 24.0, 30.0, 9.0},
		{ 4.0, 10.0, 12.0, 7.0},
		{14.0, 16.0, 36.0, 1.0},
	};
	multi::array<T, 1> const x = {1.1, 2.1, 3.1, 4.1};  // NOLINT(readability-identifier-length) BLAS naming
	{
		multi::array<T, 1> y(multi::extensions_t<1>{multi::iextension{size(a)}});  // NOLINT(readability-identifier-length) BLAS naming
		blas::gemv_n(1.0, begin(a), size(a), begin(x), 0.0, begin(y));
		BOOST_REQUIRE_CLOSE(y[1], 91.3, 0.0001);
		if(!std::is_same_v<T, float>) {  // workaround Apple Accelerate BLAS bug in dot
			BOOST_REQUIRE_CLOSE(y[2], +blas::dot(a[2], x), 0.0001);
		}
	}
	{
		multi::array<T, 1>       y(multi::extensions_t<1>{multi::iextension{size(a)}});  // NOLINT(readability-identifier-length) BLAS naming
		multi::array<T, 2> const aT = ~a;
		blas::gemv_n(1.0, begin(~aT), size(~aT), begin(x), 0.0, begin(y));
		BOOST_REQUIRE_CLOSE(y[1], 91.3, 0.0001);

		if(!std::is_same_v<T, float>) {  // workaround Apple Accelerate BLAS bug in dot
			BOOST_REQUIRE_CLOSE(y[2], +blas::dot(a[2], x), 0.0001);
		}
	}
	{
		multi::array<T, 1> y(multi::extensions_t<1>{multi::iextension{size(a)}});  // NOLINT(readability-identifier-length) BLAS naming
		auto               mv = blas::gemv(1.0, a, x);
		copy_n(mv.begin(), mv.size(), y.begin());
		BOOST_REQUIRE_CLOSE(y[1], 91.3, 0.00001);

		multi::array<T, 1> w2(multi::extensions_t<1>{multi::iextension{size(a)}});
		MV(a, x, w2);
		BOOST_REQUIRE_CLOSE(w2[0], y[0], 0.00001);
	}
	{
		multi::array<T, 1> y(multi::extensions_t<1>{multi::iextension{size(a)}});  // NOLINT(readability-identifier-length) BLAS naming
		y = blas::gemv(1.0, a, x);
		BOOST_REQUIRE_CLOSE(y[1], 91.3, 0.00001);
	}
	{
		multi::array<T, 1> y = blas::gemv(1.0, a, x);  // NOLINT(readability-identifier-length) BLAS naming
		BOOST_REQUIRE_CLOSE(y[1], 91.3, 0.00001);
	}
	{
		multi::array<T, 1> y(multi::extensions_t<1>{multi::iextension{size(a)}}, 0.);  // NOLINT(readability-identifier-length) BLAS naming
		y += blas::gemv(1.0, a, x);
		BOOST_REQUIRE_CLOSE(y[1], 91.3, 0.00001);
	}
	{
		multi::array<T, 1> y = {4.0, 5.0, 6.0};  // NOLINT(readability-identifier-length) BLAS naming
		blas::gemv(1.1, a, x, 1.0, y);  // y = a*M*x + b*y
		BOOST_REQUIRE_CLOSE(y[1], 105.43, 0.00001);
	}
}

BOOST_AUTO_TEST_CASE_TEMPLATE(multi_blas_gemv_real, T, fp_types) {
	namespace blas = multi::blas;

	using std::abs;
	// NOLINTNEXTLINE(readability-identifier-length) BLAS naming
	multi::array<T, 2> const a = {
		{ 9.0, 24.0, 30.0, 9.0},
		{ 4.0, 10.0, 12.0, 7.0},
		{14.0, 16.0, 36.0, 1.0},
	};
	multi::array<T, 1> const x = {1.1, 2.1, 3.1, 4.1};  // NOLINT(readability-identifier-length) BLAS naming
	{
		multi::array<T, 1> y     = {4.0, 5.0, 6.0};  // NOLINT(readability-identifier-length) BLAS naming
		T const            alpha = 1.1;
		T const            beta  = 1.2;
		blas::gemv(alpha, a, x, beta, y);  // y = a*M*x + b*y

		multi::array<T, 1> const y3 = {214.02, 106.43, 188.37};
		BOOST_REQUIRE( abs(y[1] - y3[1]) < 2e-14 );
	}
	if constexpr(!std::is_same_v<T, float>) {
		auto Y = +blas::gemv(1.0, a, x);  // NOLINT(readability-identifier-length) BLAS naming
		BOOST_REQUIRE_CLOSE(Y[0], +blas::dot(a[0], x), 0.00001);
		BOOST_REQUIRE_CLOSE(Y[1], +blas::dot(a[1], x), 0.00001);
		BOOST_REQUIRE_CLOSE(Y[2], +blas::dot(a[2], x), 0.00001);
	}
	{
		multi::array<T, 1> const x   = {1.0, 2.0, 3.0};  // NOLINT(readability-identifier-length) BLAS naming
		multi::array<T, 1> const y   = {4.0, 5.0, 6.0};  // NOLINT(readability-identifier-length) BLAS naming
		multi::array<T, 1> const dot = blas::gemv(1., multi::array<T, 2>({x}), y);
		if(!std::is_same_v<T, float>) {  // workaround Apple Accelerate BLAS bug in dot
			BOOST_REQUIRE( dot[0] == blas::dot(x, y) );
		}
	}
	{
		using blas::operators::operator%;
		using blas::operators::operator-;
		using blas::operators::operator^;
		BOOST_REQUIRE_SMALL(((~+~a) % x - a % x) ^ 2, 1e-9);
	}
}

BOOST_AUTO_TEST_CASE_TEMPLATE(multi_blas_gemv_real_complex, T, fp_types) {
	namespace blas = multi::blas;
	using complex  = std::complex<T>;
	using std::abs;

	// NOLINTNEXTLINE(readability-identifier-length) BLAS naming
	multi::array<complex, 2> const M = {
		{ {9.0, 0.0}, {24.0, 0.0}, {30.0, 0.0}, {9.0, 0.0}},
		{ {4.0, 0.0}, {10.0, 0.0}, {12.0, 0.0}, {7.0, 0.0}},
		{{14.0, 0.0}, {16.0, 0.0}, {36.0, 0.0}, {1.0, 0.0}},
	};

	// NOLINTNEXTLINE(readability-identifier-length) BLAS naming
	multi::array<complex, 1> const X = {
		{1.1, 0.0},
		{2.1, 0.0},
		{3.1, 0.0},
		{4.1, 0.0},
	};
	{
		// NOLINTNEXTLINE(readability-identifier-length) BLAS naming
		multi::array<complex, 1> Y = {
			{4.0, 0.0},
			{5.0, 0.0},
			{6.0, 0.0},
		};

		auto const alpha = T{1.1};
		auto const beta  = T{1.2};

		blas::gemv(alpha, M, X, beta, Y);  // y = a*M*x + b*y

		multi::array<complex, 1> const Y3 = {
			{214.02, 0.0},
			{106.43, 0.0},
			{188.37, 0.0},
		};

		using blas::operators::operator-;
		T const n2{blas::nrm2(Y - Y3)};
		BOOST_REQUIRE_SMALL(n2, T{1.0e-4});
	}
}

BOOST_AUTO_TEST_CASE(multi_blas_gemv_complex) {
	namespace blas = multi::blas;
	using complex  = std::complex<double>;
	auto const I   = complex{0.0, 1.0};  // NOLINT(readability-identifier-length) imag unit

	using std::abs;

	// NOLINTNEXTLINE(readability-identifier-length) BLAS naming
	multi::array<complex, 2> const a = {
		{2.0 + 3.0 * I, 2.0 + 1.0 * I, 1.0 + 2.0 * I},
		{4.0 + 2.0 * I, 2.0 + 4.0 * I, 3.0 + 1.0 * I},
		{7.0 + 1.0 * I, 1.0 + 5.0 * I, 0.0 + 3.0 * I},
	};
	multi::array<complex, 1> const x = {1.0 + 2.0 * I, 2.0 + 1.0 * I, 9.0 + 2.0 * I};  // NOLINT(readability-identifier-length) BLAS naming
	BOOST_REQUIRE(( +blas::gemv(1., a, x) == multi::array<complex, 1>{4.0 + 31.*I, 25.0 + 35.0*I, -4.0 + 53.0*I} ));

	auto aT = +~a;
	BOOST_REQUIRE(( +blas::gemv(1., ~aT, x) == multi::array<complex, 1>{4.0 + 31.0*I, 25.0 + 35.0*I, -4.0 + 53.0*I} ));

	BOOST_REQUIRE( +blas::gemv(1., ~a, x) == (multi::array<complex, 1>{63.0 + 38.0*I, -1.0 + 62.0*I, -4.0 + 36.0*I}) );
	BOOST_REQUIRE( +blas::gemv(1., ~a, x) == + blas::gemv(1.0, aT, x) );
}

BOOST_AUTO_TEST_CASE(multi_blas_gemv_temporary) {
	using complex = std::complex<double>;

	// NOLINTNEXTLINE(readability-identifier-length) BLAS naming
	multi::array<complex, 2> const A = {
		{{1.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}},
		{{0.0, 0.0}, {1.0, 0.0}, {0.0, 0.0}},
		{{0.0, 0.0}, {0.0, 0.0}, {1.0, 0.0}},
	};

	auto const B = [](auto array) {  // NOLINT(readability-identifier-length) BLAS naming
		auto rand = [gauss = std::normal_distribution<>{}, gen = std::mt19937{1}]() mutable { return complex{gauss(gen), gauss(gen)}; };  // NOLINT(cert-msc32-c,cert-msc51-cpp) test purposes
		std::generate(array.elements().begin(), array.elements().end(), rand);
		return array;
	}(multi::array<complex, 2>({3, 3}));

	using blas::operators::operator*;
	using blas::operators::operator-;
	using blas::operators::operator^;
	BOOST_REQUIRE( (((+(A*B))[0] - B[0])^2) == 0.0 );
	BOOST_REQUIRE( (((+(A*B))[1] - B[1])^2) == 0.0 );
	BOOST_REQUIRE( (((+(A*B))[2] - B[2])^2) == 0.0 );
}

BOOST_AUTO_TEST_CASE(multi_blas_gemv_context) {
	// NOLINTNEXTLINE(readability-identifier-length) BLAS naming
	multi::array<double, 2> const a = {
		{ 9.0, 24.0, 30.0, 9.0},
		{ 4.0, 10.0, 12.0, 7.0},
		{14.0, 16.0, 36.0, 1.0},
	};
	multi::array<double, 1> const x = {1.1, 2.1, 3.1, 4.1};  // NOLINT(readability-identifier-length) conventional name in BLAS

	blas::context ctxt;
	{
		multi::array<double, 1> y(multi::extensions_t<1>{multi::iextension{size(a)}});  // NOLINT(readability-identifier-length) conventional name in BLAS
		blas::gemv_n(&ctxt, 1.0, begin(a), size(a), begin(x), 0.0, begin(y));
		BOOST_REQUIRE_CLOSE(y[1], 91.3, 0.0001);
		BOOST_REQUIRE_CLOSE(y[2], +blas::dot(a[2], x), 0.0001);
	}
	{
		multi::array<double, 1>       y(multi::extensions_t<1>{multi::iextension{size(a)}});  // NOLINT(readability-identifier-length) conventional name in BLAS
		multi::array<double, 2> const aT = ~a;
		blas::gemv_n(&ctxt, 1.0, begin(~aT), size(~aT), begin(x), 0.0, begin(y));
		BOOST_REQUIRE_CLOSE(y[1], 91.3, 0.00001);
		BOOST_REQUIRE_CLOSE(y[2], +blas::dot(a[2], x), 0.00001);
	}
	{
		multi::array<double, 1> y(multi::extensions_t<1>{multi::iextension{size(a)}});  // NOLINT(readability-identifier-length) conventional name in BLAS
		auto&&                  mv = blas::gemv(&ctxt, 1.0, a, x);
		copy_n(mv.begin(), mv.size(), y.begin());
		BOOST_REQUIRE_CLOSE(y[1], 91.3, 0.00001);
	}
	{
		multi::array<double, 1> y(multi::extensions_t<1>{multi::iextension{size(a)}});  // NOLINT(readability-identifier-length) conventional name in BLAS
		y = blas::gemv(&ctxt, 1.0, a, x);
		BOOST_REQUIRE_CLOSE(y[1], 91.3, 0.00001);
	}
	{
		multi::array<double, 1> y(multi::extensions_t<1>{multi::iextension{size(a)}});  // NOLINT(readability-identifier-length) conventional name in BLAS
		y = blas::gemv(1.0, a, x);
		BOOST_REQUIRE_CLOSE(y[1], 91.3, 0.00001);
	}
	{
		multi::array<double, 1> y(multi::extensions_t<1>{multi::iextension{size(a)}});  // NOLINT(readability-identifier-length) conventional name in BLAS
		y() = blas::gemv(1.0, a, x);
		BOOST_REQUIRE_CLOSE(y[1], 91.3, 0.00001);
	}
	{
		multi::array<double, 1> y = blas::gemv(&ctxt, 1.0, a, x);  // NOLINT(readability-identifier-length) conventional name in BLAS
		BOOST_REQUIRE_CLOSE(y[1], 91.3, 0.00001);
	}
	{
		multi::array<double, 1> y(multi::extensions_t<1>{multi::iextension{size(a)}}, 0.0);  // NOLINT(readability-identifier-length) conventional name in BLAS
		y += blas::gemv(&ctxt, 1.0, a, x);
		BOOST_REQUIRE_CLOSE(y[1], 91.3, 0.00001);
	}
	{
		multi::array<double, 1> y = {4.0, 5.0, 6.0};  // NOLINT(readability-identifier-length) conventional name in BLAS
		y += blas::gemv(&ctxt, 1.1, a, x);
		BOOST_REQUIRE_CLOSE(y[1], 105.43, 0.00001);
	}
}
