// Copyright 2025 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/multi/array.hpp>  // for array, implicit_cast, explicit_cast

#include <boost/core/lightweight_test.hpp>

#include <algorithm>  // for copy
#include <cmath>
#include <iostream>
#include <iterator>  // for distance
#include <vector>

namespace multi = boost::multi;

class uniform_cspline {  // NOLINT(misc-use-internal-linkage)
	using argument_type = double;
	using result_type   = double;  // typename std::decay<decltype(std::declval<vector>()[0])>::type;
	using vector        = std::vector<result_type>;
	using size_type     = vector::size_type;
	using index         = multi::array<double, 2>::index;

	argument_type lower_;
	argument_type dx_;

 public:
	auto dx() const -> argument_type { return dx_; }
	auto lower() const -> argument_type { return lower_; }
	auto upper() const -> argument_type { return lower_ + (dx_ * static_cast<double>(K_.size() - 1)); }

 private:
	multi::array<double, 2> K_;
	// vector a;
	// vector b;
	// vector c;
	// vector d;
	// std::vector<std::array<double, 4>> K;
 public:
	// template<class It>
	// uniform_cspline(It a_begin, It a_end, argument_type lower, argument_type dx)
	// : uniform_cspline(std::vector<double>(a_begin, a_end), lower, dx) {}

	template<class It>
	uniform_cspline(It a_first, It a_last, argument_type lower, argument_type dx)  // NOLINT(bugprone-easily-swappable-parameters)
	: lower_{lower}, dx_{dx}, K_({static_cast<multi::array<double, 2>::size_type>(std::distance(a_first, a_last)), 4}) {
		auto const n = K_.size();        // NOLINT(readability-identifier-length)
		auto&&     a = K_.rotated()[0];  // NOLINT(readability-identifier-length)

		std::copy(a_first, a_last, a.begin());

		auto&& b = K_.rotated()[1];  // NOLINT(readability-identifier-length)
		auto&& c = K_.rotated()[2];  // NOLINT(readability-identifier-length)
		auto&& d = K_.rotated()[3];  // NOLINT(readability-identifier-length)

		auto&& A = b;  // NOLINT(readability-identifier-length)
		auto&& l = c;  // NOLINT(readability-identifier-length)
		auto&& u = d;  // NOLINT(readability-identifier-length)
		auto&& z = A;  // NOLINT(readability-identifier-length)

		for(index i = 1; i != n - 1; ++i) {  // NOLINT(altera-unroll-loops,altera-id-dependent-backward-branch) TODO(correaa) use algorithms
			A[i] = (3 * (a[i + 1] - a[i]) / dx_) - (3 * (a[i] - a[i - 1]) / dx_);
		}

		l[0] = 4 * dx_;
		u[0] = dx_ / l[0];
		z[0] = A[1] / l[0];

		for(index i = 1; i != n - 2; ++i) {  // NOLINT(altera-unroll-loops,altera-id-dependent-backward-branch) TODO(correaa) use algorithms
			l[i] = (4 * dx_) - (dx_ * u[i - 1]);
			u[i] = dx_ / l[i];
			z[i] = (A[i + 1] - (dx_ * z[i - 1])) / l[i];
		}

		c[n - 1] = 0;

		for(index j = n - 2; j != 0; --j) {  // NOLINT(altera-unroll-loops,altera-id-dependent-backward-branch) TODO(correaa) use algorithms
			c[j] = z[j - 1] - (u[j - 1] * c[j + 1]);
			b[j] = ((a[j + 1] - a[j]) / dx_) - (dx_ * (c[j + 1] + (2 * c[j])) / 3);
			d[j] = (c[j + 1] - c[j]) / (3 * dx_);
		}

		c[0] = 0;
		b[0] = ((a[1] - a[0]) / dx_) - (dx_ * c[1] / 3);
		d[0] = c[1] / (3 * dx_);
		// for(std::size_t i = 0; i != K.size(); ++i) K[i] = {a[i], b[i], c[i], d[i]};
	}
	auto operator()(argument_type x) const -> result_type {
		auto const i   = static_cast<index>((x - lower_) / dx_);
		auto const Dx  = x - (static_cast<double>(i) * dx_) - lower_;
		// return a[i] + Dx*(b[i] + Dx*(c[i] + Dx*d[i]));
		// auto const& Ki = K[i]; using std::get;
		// return K[0][i] + Dx*(K[1][i] + Dx*(K[2][i] + Dx*K[3][i]));
		auto const& Ki = K_[i];
		return Ki[0] + (Dx * (Ki[1] + (Dx * (Ki[2] + (Dx * Ki[3])))));
	}
};

// #if defined(__cpp_deduction_guides)
// template<class It, class ValueType> uniform_cspline(It, It, ValueType, ValueType) -> uniform_cspline<typename std::iterator_traits<It>::value_type>;
// #endif

auto main() -> int {  // NOLINT(readability-function-cognitive-complexity,bugprone-exception-escape)
	using std::cout;
	using vector = std::vector<double>;

	//           0.0,  0.4, 0.8, 1.2,  1.6
	vector a = {10.0, 12.0, 9.0, 8.0, 13.0};  // NOLINT(readability-identifier-length)

	uniform_cspline const as(a.begin(), a.end(), 1.0, 0.4);
	{
		// for(double x = as.lower() + 0.01; x < as.upper(); x += 0.01) {
		for(int i = 0; i != static_cast<int>((as.upper() - as.lower()) / 0.01); ++i) {  // NOLINT(altera-unroll-loops,altera-id-dependent-backward-branch)
			double const x = as.lower() + (static_cast<double>(i) * 0.01);              // NOLINT(readability-identifier-length)
			std::cout << x << '\t' << as(x) << '\n';
		}
	}

	BOOST_TEST( std::abs(as.dx() - 0.4) < 1.0e-10 );
	BOOST_TEST( std::abs( as.lower() - 1.0) < 1.0e-10 );
	BOOST_TEST( std::abs( as.upper() - ( 1.0 + (as.dx()*static_cast<double>(a.size() - 1))) ) < 1.0e-10 );
	BOOST_TEST( std::abs( as(1.0) - 10.0) < 1.0e-10 );
	BOOST_TEST( std::abs( as(1.4) - 12.0) < 1.0e-10 );

	{
		// using Clock  = std::chrono::high_resolution_clock;
		// using ns     = std::chrono::nanoseconds;

		// std::mt19937 eng{std::random_device{}()};
		// std::uniform_real_distribution<double>{-10.0, 10.0};

		// auto gen = [
		// 	dist = std::uniform_real_distribution<>{-10.0, 10.0},
		// 	eng = std::mt19937{std::random_device{}()}
		// ]() mutable { return dist(eng); };

		// int const N = 1000;

		// std::vector<vector> rep(100000);

		// generate(rep.begin(), rep.end(), [&]() { vector ret(N); generate(ret.begin(), ret.end(), gen); return ret; });
		// {
		// 	cout << "using val semantics\nconstruction: repeats: " << rep.size() << " size: " << rep[0].size() << '\n';
		// 	auto tic = Clock::now();
		// 	for(std::size_t i = 0; i != rep.size(); ++i) {
		// 		uniform_cspline ucsp(rep[i].begin(), rep[i].end(), 1.0, 0.4);
		// 	}
		// 	auto toc = Clock::now();
		// 	cout << '\t' << static_cast<double>(ns{toc - tic}.count()) / 1.0e9 << " sec\n";
		// }
		// {
		// 	vector a(10000);
		// 	generate(begin(a), end(a), gen);
		// 	uniform_cspline as(begin(a), end(a), 0., 0.4);
		// 	{
		// 		std::uniform_real_distribution<double> x_dist{as.left_endpoint(), as.right_endpoint()};
		// 		vector                                 x(1000000000);
		// 		double                                 nodiscard = 0;
		// 		generate(begin(x), end(x), [&x_dist, &eng]() { return x_dist(eng); });
		// 		{
		// 			cout << "evaluation: repeats: " << x.size() << " for size: " << a.size() << '\n';
		// 			auto tic = Clock::now();
		// 			for(std::size_t i = 0; i != x.size(); ++i)
		// 				nodiscard += as(x[i]);
		// 			auto toc = Clock::now();
		// 			cout << '\t' << ns{toc - tic}.count() / 1e9 << " sec [" << (char)nodiscard << "]\n";
		// 		}
		// 	}
		// }
	}

	return boost::report_errors();
}
