// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// © Alfredo Correa 2021

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi gemm (not blas)"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include "../array.hpp"

#include "../algorithms/gemm.hpp"

#include <numeric>
#include <random>

namespace multi = boost::multi;

BOOST_AUTO_TEST_CASE(algorithm_gemm) {

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-1.0, +1.0);

	auto const N = 256;

	auto const A = [&] {
		multi::array<double, 2> _({N, N});
		std::generate(begin(elements(_)), end(elements(_)), [&]{return dis(gen);});
		return _;
	}();

	auto const B = [&] {
		multi::array<double, 2> _({(~A).size(), N});
		std::generate(begin(elements(_)), end(elements(_)), [&]{return dis(gen);});
		return _;
	}();

	// zero init, beta = zero multiplication
	{
		multi::array<double, 2> C_gold({  A.size() , (~B).size()}, 0.);
		multi::array<double, 2> C = C_gold;

		multi::detail::naive_gemm(1., A, B, 0., C_gold);
		multi::gemm(1., A, B, 0., C);

		BOOST_TEST( C[123][121] == C_gold[123][121] , boost::test_tools::tolerance(1e-12) );
	}

	// non-zero init, beta = zero multiplication
	{
		multi::array<double, 2> C_gold({  A.size() , (~B).size()}, 0.);
		std::generate(begin(elements(C_gold)), end(elements(C_gold)), [&]{return dis(gen);});

		multi::array<double, 2> C = C_gold;

		multi::detail::naive_gemm(1., A, B, 0., C_gold);
		multi::gemm(1., A, B, 0., C);

		BOOST_TEST( C[123][121] == C_gold[123][121] , boost::test_tools::tolerance(1e-12) );
	}
	// non-zero init, beta = one multiplication
	{
		multi::array<double, 2> C_gold({  A.size() , (~B).size()}, 0.);
		std::generate(begin(elements(C_gold)), end(elements(C_gold)), [&]{return dis(gen);});

		multi::array<double, 2> C = C_gold;

		multi::detail::naive_gemm(1., A, B, 1., C_gold);
		multi::gemm(1., A, B, 1., C);

		BOOST_TEST( C[123][121] == C_gold[123][121] , boost::test_tools::tolerance(1e-12) );
	}
	{
		multi::array<double, 2> C_gold({  A.size() , (~B).size()}, 0.);
		std::generate(begin(elements(C_gold)), end(elements(C_gold)), [&]{return dis(gen);});

		multi::array<double, 2> C = C_gold;

		multi::detail::naive_gemm(1., A, B, 0.3, C_gold);
		multi::gemm(1., A, B, 0.3, C);

		BOOST_TEST( C[123][121] == C_gold[123][121] , boost::test_tools::tolerance(1e-12) );
	}
	{
		multi::array<double, 2> C_gold({  A.size() , (~B).size()}, 0.);
		std::generate(begin(elements(C_gold)), end(elements(C_gold)), [&]{return dis(gen);});

		multi::array<double, 2> C = C_gold;

		multi::detail::naive_gemm(2., A, B, 0., C_gold);
		multi::gemm(2., A, B, 0., C);

		BOOST_TEST( C[123][121] == C_gold[123][121] , boost::test_tools::tolerance(1e-12) );
	}
	{
		multi::array<double, 2> C_gold({  A.size() , (~B).size()}, 0.);
		std::generate(begin(elements(C_gold)), end(elements(C_gold)), [&]{return dis(gen);});

		multi::array<double, 2> C = C_gold;

		multi::detail::naive_gemm(2., A, B, 0.3, C_gold);
		multi::gemm(2., A, B, 0.3, C);

		BOOST_TEST( C[123][121] == C_gold[123][121] , boost::test_tools::tolerance(1e-12) );
	}
}
