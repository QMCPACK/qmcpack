// Copyright 2019-2025 Alfredo A. Correa

#include <boost/core/lightweight_test.hpp>

#include <boost/multi/adaptors/lapack/core.hpp>
#include <boost/multi/adaptors/lapack/getrf.hpp>

#include <boost/multi/array.hpp>

#include <boost/multi/adaptors/blas/gemm.hpp>
#include <boost/multi/adaptors/blas/gemv.hpp>

#include <array>
#include <cassert>
#include <cmath>
#include <iostream>

namespace multi = boost::multi;

// BOOST_AUTO_TEST_CASE(lapack_getrf){

//// https://www.ibm.com/support/knowledgecenter/SSFHY8_6.2/reference/am5gr_hsgetrf.html
//  multi::array<double, 2> A = {
//      { 1.0,  1.2,  1.4,  1.6,  1.8,  2.0,  2.2,  2.4,  2.6 },
//      { 1.2,  1.0,  1.2,  1.4,  1.6,  1.8,  2.0,  2.2,  2.4 },
//      { 1.4,  1.2,  1.0,  1.2,  1.4,  1.6,  1.8,  2.0,  2.2 },
//      { 1.6,  1.4,  1.2,  1.0,  1.2,  1.4,  1.6,  1.8,  2.0 },
//      { 1.8,  1.6,  1.4,  1.2,  1.0,  1.2,  1.4,  1.6,  1.8 },
//      { 2.0,  1.8,  1.6,  1.4,  1.2,  1.0,  1.2,  1.4,  1.6 },
//      { 2.2,  2.0,  1.8,  1.6,  1.4,  1.2,  1.0,  1.2,  1.4 },
//      { 2.4,  2.2,  2.0,  1.8,  1.6,  1.4,  1.2,  1.0,  1.2 },
//      { 2.6,  2.4,  2.2,  2.0,  1.8,  1.6,  1.4,  1.2,  1.0 }
//  };

//  multi::array<int, 1> P({9}, 0.);
//  lapack::context ctxt;
//  auto const& LU = multi::lapack::getrf(ctxt, A, P);

//  BOOST_REQUIRE( LU.size() == A.size() );

//  BOOST_REQUIRE_CLOSE( LU[0][0] , 2.6      , 1e-5 );
//  BOOST_REQUIRE_CLOSE( LU[0][8] , 1.       , 1e-5 );
//  BOOST_REQUIRE_CLOSE( LU[8][0] , 0.923077 , 1e-5 );
//  BOOST_REQUIRE_CLOSE( LU[8][8] , 0.4      , 1e-5 );

//}

// BOOST_AUTO_TEST_CASE(lapack_getrf2){

//// https://www.ibm.com/support/knowledgecenter/SSFHY8_6.2/reference/am5gr_hsgetrf.html
//  multi::array<double, 2> A = {
//      { 1.0,  1.0,  1.0,  1.0,  0.0,  0.0,   0.0,   0.0,   0.0 },
//      { 1.0,  1.0,  1.0,  1.0,  1.0,  0.0,   0.0,   0.0,   0.0 },
//      { 4.0,  1.0,  1.0,  1.0,  1.0,  1.0,   0.0,   0.0,   0.0 },
//      { 0.0,  5.0,  1.0,  1.0,  1.0,  1.0,   1.0,   0.0,   0.0 },
//      { 0.0,  0.0,  6.0,  1.0,  1.0,  1.0,   1.0,   1.0,   0.0 },
//      { 0.0,  0.0,  0.0,  7.0,  1.0,  1.0,   1.0,   1.0,   1.0 },
//      { 0.0,  0.0,  0.0,  0.0,  8.0,  1.0,   1.0,   1.0,   1.0 },
//      { 0.0,  0.0,  0.0,  0.0,  0.0,  9.0,   1.0,   1.0,   1.0 },
//      { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  10.0,  11.0,  12.0 }
//  };

//  multi::array<int, 1> P({9}, 0.);
//  lapack::context ctxt;
//  auto const& LU = multi::lapack::getrf(ctxt, A, P);

//  BOOST_REQUIRE( LU.size() == A.size() );

//  for(int i = 0; i != 9; ++i){
//      for(int j = 0; j != 9; ++j){
//          std::cout<<'\t'<< LU[i][j] <<',';
//      }
//      std::cout<<std::endl;
//  }
//  std::cout<<std::endl;
//  for(int i = 0; i != 9; ++i){std::cout<< P[i] <<", ";}
//  std::cout<<std::endl;

//}

auto main() -> int {  // NOLINT(bugprone-exception-escape)
	// BOOST_AUTO_TEST_CASE(lapack_getrf)
	{
		multi::array<double, 2> const Aconst = {
			{ 6.80, -6.05, -0.45,  8.32, -9.67},
			{-2.11, -3.30,  2.58,  2.71, -5.14},
			{ 5.66,  5.36, -2.70,  4.35, -7.26},
			{ 5.97, -4.44,  0.27, -7.17,  6.08},
			{ 8.23,  1.08,  9.04,  2.14, -6.87}
		};

		multi::array<double, 2> const Bconst = {
			{ 4.02, -1.56,  9.81},
			{ 6.19,  4.00, -4.09},
			{-8.22, -8.67, -4.57},
			{-7.57,  1.75, -8.61},
			{-3.03,  2.86,  8.99}
		};

		multi::array<multi::lapack::index, 1> PP({5}, 0.0);

		auto AA = Aconst;
		auto BB = Bconst;

		auto AT = +~AA;
		auto BT = +~BB;

		auto lu_solve = [](auto&& Aio, auto&& Po, auto&& Bio) {  // solve A.X = B; put result in B
			auto AT = +~Aio;
			auto BT = +~Bio;

			auto const& LU = multi::lapack::getrf(~AT, Po);
			assert(LU.size() == Po.size());
			multi::lapack::getrs(LU, std::as_const(Po), ~BT);

			Bio = ~BT;
			Aio = ~AT;
		};
		lu_solve(~AT, PP, ~BT);

		// using multi::blas::operators::operator*;
		BOOST_TEST( std::abs((+multi::blas::gemm(1.0, Aconst, (~BT)))[1][2] - Bconst[1][2]) < 1e-10);
	}

	// BOOST_AUTO_TEST_CASE(lapack_getrf_two_column)
	{
		multi::array<double, 2> const Aconst = {
			{ 6.80, -6.05, -0.45,  8.32, -9.67},
			{-2.11, -3.30,  2.58,  2.71, -5.14},
			{ 5.66,  5.36, -2.70,  4.35, -7.26},
			{ 5.97, -4.44,  0.27, -7.17,  6.08},
			{ 8.23,  1.08,  9.04,  2.14, -6.87}
		};

		multi::array<double, 2> const Bconst = {
			{ 4.02, -1.56},
			{ 6.19,  4.00},
			{-8.22, -8.67},
			{-7.57,  1.75},
			{-3.03,  2.86}
		};

		multi::array<multi::lapack::index, 1> PP({5}, 0.0);

		auto AA = Aconst;
		auto BB = Bconst;

		auto AT = +~AA;
		auto BT = +~BB;

		auto lu_solve = [](auto&& Aio, auto&& Po, auto&& Bio) {  // solve A.X = B; put result in B
			auto AT = +~Aio;
			auto BT = +~Bio;

			auto const& LU = multi::lapack::getrf(~AT, Po);
			assert(LU.size() == Po.size());
			multi::lapack::getrs(LU, std::as_const(Po), ~BT);

			Bio = ~BT;
			Aio = ~AT;
		};
		lu_solve(~AT, PP, ~BT);

		// using multi::blas::operators::operator*;
		// BOOST_REQUIRE_CLOSE( (Aconst*(~BT))[2][1] , Bconst[2][1] , 1e-10);
		// BOOST_REQUIRE_CLOSE( (Aconst*(~BT))[2][0] , Bconst[2][0] , 1e-10);
	}

	// BOOST_AUTO_TEST_CASE(lapack_getrf_one_column)
	{
		multi::array<double, 2> const Aconst = {
			{ 6.80, -6.05, -0.45,  8.32, -9.67},
			{-2.11, -3.30,  2.58,  2.71, -5.14},
			{ 5.66,  5.36, -2.70,  4.35, -7.26},
			{ 5.97, -4.44,  0.27, -7.17,  6.08},
			{ 8.23,  1.08,  9.04,  2.14, -6.87}
		};

		multi::array<double, 2> const Bconst = {
			{4.02},
			{6.19},
			{-8.22},
			{-7.57},
			{-3.03}
		};

		multi::array<multi::lapack::index, 1> PP({5}, 0.0);

		auto AA = Aconst;
		auto BB = Bconst;

		auto lu_solve = [](auto&& Aio, auto&& Po, auto&& Bio) {  // solve A.X = B; put result in B
			auto AT = +~Aio;
			auto BT = +~Bio;

			auto const& LU = multi::lapack::getrf(~AT, Po);
			assert(LU.size() == Po.size());
			multi::lapack::getrs(LU, std::as_const(Po), ~BT);

			Bio = ~BT;
			Aio = ~AT;
		};
		lu_solve(AA, PP, BB);

		// using multi::blas::operators::operator*;
		// BOOST_REQUIRE_CLOSE( (Aconst*B)[1][0] , Bconst[1][0] , 1e-10);
		// BOOST_REQUIRE_CLOSE( (Aconst*B)[1][0] , Bconst[1][0] , 1e-10);
	}

	// BOOST_AUTO_TEST_CASE(lapack_getrf_one_vector)
	{
		multi::array<double, 2> const Aconst = {
			{ 6.80, -6.05, -0.45,  8.32, -9.67},
			{-2.11, -3.30,  2.58,  2.71, -5.14},
			{ 5.66,  5.36, -2.70,  4.35, -7.26},
			{ 5.97, -4.44,  0.27, -7.17,  6.08},
			{ 8.23,  1.08,  9.04,  2.14, -6.87}
		};

		multi::array<double, 1> Vconst = {4.02, 6.19, -8.22, -7.57, -3.03};

		multi::array<multi::lapack::index, 1> PP({5}, 0.0);

		auto AA = Aconst;
		auto VV = Vconst;

		auto lu_solve_one = [](auto&& Aio, auto&& Po, auto&& VV) {  // solve A.X = B; put result in B
			auto AT = +~Aio;
			//  auto BT = +~Bio;

			auto const& LU = multi::lapack::getrf(~AT, Po);
			assert(LU.size() == Po.size());
			multi::lapack::getrs_one(LU, std::as_const(Po), VV);

			//  Bio = ~BT;
			Aio = ~AT;
		};
		lu_solve_one(AA, PP, VV);

		using multi::blas::operators::operator%;
		BOOST_TEST( std::abs((Aconst % VV)[2] - Vconst[2]) < 1e-10 );
	}

	// BOOST_AUTO_TEST_CASE(lapack_getrs)
	{
		// https://www.ibm.com/support/knowledgecenter/SSFHY8_6.2/reference/am5gr_hsgetrf.html
		multi::array<double, 2> AA = {
			{1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6},
			{1.2, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4},
			{1.4, 1.2, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2},
			{1.6, 1.4, 1.2, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0},
			{1.8, 1.6, 1.4, 1.2, 1.0, 1.2, 1.4, 1.6, 1.8},
			{2.0, 1.8, 1.6, 1.4, 1.2, 1.0, 1.2, 1.4, 1.6},
			{2.2, 2.0, 1.8, 1.6, 1.4, 1.2, 1.0, 1.2, 1.4},
			{2.4, 2.2, 2.0, 1.8, 1.6, 1.4, 1.2, 1.0, 1.2},
			{2.6, 2.4, 2.2, 2.0, 1.8, 1.6, 1.4, 1.2, 1.0}
		};

		multi::array<int, 1> PP({9}, 0.0);

		auto BB = +~multi::array<double, 2>{
			{93.0, 186.0},  //  279.0,  372.0,  465.0 },
			{84.4, 168.8},  //  253.2,  337.6,  422.0 },
			{76.6, 153.2},  //  229.8,  306.4,  383.0 },
			{70.0, 140.0},  //  210.0,  280.0,  350.0 },
			{65.0, 130.0},  //  195.0,  260.0,  325.0 },
			{62.0, 124.0},  //,  186.0,  248.0,  310.0 },
			{61.4, 122.8},  //  184.2,  245.6,  307.0 },
			{63.6, 127.2},  //  190.8,  254.4,  318.0 },
			{69.0, 138.0}   //,  207.0,  276.0,  345.0 }
		};

		lapack::context const ctxt;
		multi::lapack::getrf(ctxt, ~AA, PP);

		multi::array<int, 1> dee({9}, 0.0);
		for(int i = 0; i != 9; ++i) {  // NOLINT(altera-unroll-loops)
			dee[PP[i]] = i;
		}

		for(int i = 0; i != size(AA); ++i) {
			for(int j = 0; j != size(~AA); ++j) {  // NOLINT(altera-unroll-loops)
				std::cout << '\t' << AA[i][j] << ',';
			}
			std::cout << '\n';
		}

		for(int i = 0; i != size(BB); ++i) {
			for(int j = 0; j != size(~BB); ++j) {  // NOLINT(altera-unroll-loops)
				std::cout << '\t' << BB[i][j] << ',';
			}
			std::cout << '\n';
		}
	}
}
