// Copyright 2019-2025 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/multi/adaptors/lapack/geqrf.hpp>

// #include <boost/multi/adaptors/blas/gemm.hpp>
#include <boost/multi/array.hpp>

#include <boost/core/lightweight_test.hpp>

#include <algorithm>
// #include <cmath>  // for std::abs
// IWYU pragma: no_include <memory>  // for std::allocator

namespace multi = boost::multi;

// template<class Array> auto print(Array const& arr) -> decltype(auto) {
//  using std::cout;
//  using multi::size;
//  for(int i = 0; i != arr.size(); ++i) {
//      for(int j = 0; j != arr[i].size(); ++j) { cout << arr[i][j] << ' '; }  // NOLINT(altera-unroll-loops)
//      cout << '\n';
//  }
//  return cout << '\n';
// }

// template<class Array1D> auto print_1d(Array1D const& arr) -> decltype(auto) {
//  using std::cout;
//  using multi::size;
//  for(int i = 0; i != size(arr); ++i) { cout<< arr[i] <<' '; }  // NOLINT(altera-unroll-loops)
//  return cout << '\n';
// }

auto main() -> int {  // NOLINT(bugprone-exception-escape)
	multi::array<double, 2> AA = {
		{1.0, 2.0, 3.0},
		{4.0, 5.0, 6.0},
		{7.0, 8.0, 9.0}
	};
	// lapack::context ctxt;

	multi::array<double, 1> TAU(std::min(AA.size(), (~AA).size()));
	// multi::array<double, 1> WORK(std::max(1l, 3*size(A)-1));

	multi::lapack::geqrf(AA, TAU);

	// print(A);
	// print(TAU);

	// {
	//  multi::array<double, 2> const AA = {
	//      {0.5, 1.0},
	//      {2.0, 2.5},
	//  };

	//  auto const [UU, ss, VV] = multi::lapack::gesvd(AA);  // AA == UU.Diag(ss).(VV^T)

	//  multi::array<double, 2> SS({ss.size(), ss.size()}, 0.0);
	//  SS.diagonal() = ss;

	//  auto const AA_test = +multi::blas::gemm(1.0, UU, +multi::blas::gemm(1.0, SS, ~VV));  // A_test <- UU * SS * VV^T

	//  BOOST_TEST( std::abs(AA_test[0][0] - AA[0][0]) < 1.0e-4 );
	//  BOOST_TEST( std::abs(AA_test[0][1] - AA[0][1]) < 1.0e-4 );
	//  BOOST_TEST( std::abs(AA_test[1][0] - AA[1][0]) < 1.0e-4 );
	//  BOOST_TEST( std::abs(AA_test[1][1] - AA[1][1]) < 1.0e-4 );
	// }

	// {
	//  multi::array<double, 2> AA = {
	//      { 2.27,  0.94, 1.07,  0.63, -2.35,  0.62},
	//      {-1.54, -0.78, 1.22,  2.93,  2.30, -7.39},
	//      { 1.15, -0.48, 0.79, -1.45,  1.03,  1.03},
	//      {-1.94, -3.09, 0.63,  2.30, -2.57, -2.57},
	//  };

	//  auto const AA_copy = AA;

	//  // Output arrays
	//  multi::array<double, 1> ss(std::min(AA.size(), (~AA).size()));  // Singular values

	//  multi::array<double, 2> UU({(~AA).size(), (~AA).size()});  // Left singular vectors
	//  multi::array<double, 2> VT({AA.size(), AA.size()});        // Right singular vectors

	//  boost::multi::lapack::gesvd(AA, UU, ss, VT);

	//  std::cout << "Original array:\n";
	//  {
	//      auto [is, js] = AA.extensions();
	//      for(auto i : is) {
	//          for(auto j : js) {  // NOLINT(altera-unroll-loops)
	//              std::cout << AA_copy[i][j] << ' ';
	//          }
	//          std::cout << '\n';
	//      }
	//  }

	//  multi::array<double, 2> SS({ss.extension(), ss.extension()}, 0.0);
	//  std::copy(ss.begin(), ss.end(), SS.diagonal().begin());

	//  // Print singular values
	//  std::cout << "Singular values:\n";
	//  for(auto i : ss.extension()) {  // NOLINT(altera-unroll-loops)
	//      std::cout << ss[i] << ' ';
	//  }
	//  std::cout << '\n';

	//  std::cout << "Singular vectors as array:\n";
	//  for(auto const& row : SS) {
	//      for(auto const& elem : row) {  // NOLINT(altera-unroll-loops)
	//          std::cout << elem << ' ';
	//      }
	//      std::cout << '\n';
	//  }

	//  // Print left singular vectors
	//  std::cout << "Left singular vectors:\n";
	//  for(auto const& row : UU) {
	//      for(auto const& elem : row) {  // NOLINT(altera-unroll-loops)
	//          std::cout << elem << ' ';
	//      }
	//      std::cout << '\n';
	//  }

	//  // Print right singular vectors
	//  std::cout << "Right singular vectors:\n";
	//  {
	//      auto [is, js] = VT.extensions();
	//      for(auto i : is) {
	//          for(auto j : js) {  // NOLINT(altera-unroll-loops)
	//              std::cout << VT[i][j] << ' ';
	//          }
	//          std::cout << '\n';
	//      }
	//  }
	//  // Singular values: s
	//  // 9.43397 4.71924 3.19716 1.88613
	//  // Left singular vectors: UU
	//  // -0.282618 -0.218724 0.114543 0.392407 0.114648 -0.831889
	//  // 0.0753653 0.370403 -0.0655296 -0.282965 0.866854 -0.146024
	//  // 0.697198 0.484235 0.360692 0.216255 -0.210349 -0.241494
	//  // 0.338129 -0.66227 0.546151 -0.338876 0.184255 -2.94997e-06
	//  // 0.558963 -0.354894 -0.725448 0.116975 0.0640779 -0.132465
	//  // 0.039885 -0.126202 0.167128 0.768545 0.391302 0.459099
	//  // Right singular vectors: VT
	//  // -0.133831 -0.393446 0.908492 0.0439558
	//  // 0.880508 0.372702 0.28873 0.0493377
	//  // -0.152352 0.213988 0.0235597 0.964595
	//  // 0.428467 -0.812713 -0.301202 0.255325

	//  // BOOST_TEST( std::abs(ss[0] - 9.43397 ) < 1e-4 );
	//  // BOOST_TEST( std::abs(ss[3] - 1.88613 ) < 1e-4 );

	//  // BOOST_TEST( std::abs(UU[0][0] - -0.282618) < 1e-4 );
	//  // BOOST_TEST( std::abs(UU[0][5] - -0.831889) < 1e-4 );
	//  // BOOST_TEST( std::abs(UU[5][0] -  0.039885) < 1e-4 );
	//  // BOOST_TEST( std::abs(UU[5][5] -  0.459099) < 1e-4 );

	//  // BOOST_TEST( std::abs(VT[0][0] - -0.133831) < 1e-4 );
	//  // BOOST_TEST( std::abs(VT[0][3] - 0.0439558) < 1e-4 );
	//  // BOOST_TEST( std::abs(VT[3][0] -  0.42846 ) < 1e-4 );
	//  // BOOST_TEST( std::abs(VT[3][3] -  0.255325) < 1e-4 );
	// }

	// {
	//  multi::array<double, 2> AA = {
	//      {0.5, 1.0},
	//      {2.0, 2.5},
	//  };

	//  auto const AA_gold = AA;

	//  // Output arrays
	//  multi::array<double, 1> ss(std::min(AA.size(), (~AA).size()));  // Singular values

	//  multi::array<double, 2> UU({(~AA).size(), (~AA).size()});  // Left singular vectors
	//  multi::array<double, 2> VT({AA.size(), AA.size()});        // Right singular vectors

	//  multi::lapack::gesvd(AA, UU, ss, VT);  // AA == VT.SS.(UU^T)

	//  multi::array<double, 2> SS({ss.extension(), ss.extension()}, 0.0);
	//  std::copy(ss.begin(), ss.end(), SS.diagonal().begin());

	//  // auto const SSUUT   = +multi::blas::gemm(1.0, SS, ~UU);
	//  auto const AA_test = +multi::blas::gemm(1.0, VT, +multi::blas::gemm(1.0, SS, ~UU));

	//  BOOST_TEST( std::abs(AA_test[0][0] - AA_gold[0][0]) < 1.0e-4 );
	//  BOOST_TEST( std::abs(AA_test[0][1] - AA_gold[0][1]) < 1.0e-4 );
	//  BOOST_TEST( std::abs(AA_test[1][0] - AA_gold[1][0]) < 1.0e-4 );
	//  BOOST_TEST( std::abs(AA_test[1][1] - AA_gold[1][1]) < 1.0e-4 );
	// }

	return boost::report_errors();
}

// BOOST_AUTO_TEST_CASE(lapack_geqrf){

//  multi::array<double, 2> A =
//      {
//          {1.0, 2.0, 3.0},
//          {4.0, 5.0, 6.0},
//          {7.0, 8.0, 9.0}
//      }
//  ;
// //  multi::lapack::context ctxt;

//  multi::array<double, 1> TAU(std::min(size(A), size(~A)));
//  multi::array<double, 1> WORK(std::max(1l, 3*size(A)-1));

//  multi::lapack::geqrf(ctxt, A, TAU, WORK);

//  print(A);
//  print(TAU);

// }

// #if 0
// BOOST_AUTO_TEST_CASE(lapack_syev, *boost::unit_test::tolerance(0.00001) ){
// {
//  multi::array<double, 2> A = {
//      {167.413, 126.804, 125.114},
//      {NAN    , 167.381, 126.746},
//      {NAN    , NAN    , 167.231}
//  };
//  multi::array<double, 1> W(size(A));
//  multi::array<double, 1> WORK(std::max(1l, 3*size(A)-1));
//  multi::lapack::syev(multi::blas::filling::upper, A, W, WORK);
//  BOOST_TEST( A[2][1] == -0.579092 );
//  BOOST_TEST( W[1] == 42.2081 );
// }
// {
//  multi::array<double, 2> A = {
//      {167.413, 126.804, 125.114},
//      {NAN    , 167.381, 126.746},
//      {NAN    , NAN    , 167.231}
//  };
//  multi::array<double, 1> W(size(A));
//  multi::lapack::syev(multi::blas::filling::upper, A, W);
//  BOOST_TEST( A[2][1] == -0.579092 );
//  BOOST_TEST( W[1] == 42.2081 );
// }
// {
//  multi::array<double, 2> A = {
//      {167.413, 126.804, 125.114},
//      {NAN    , 167.381, 126.746},
//      {NAN    , NAN    , 167.231}
//  };
//  multi::array<double, 1> W(size(A));
//  multi::lapack::syev(multi::blas::filling::lower, rotated(A), W);
//  BOOST_TEST( A[2][1] == -0.579092 );
//  BOOST_TEST( W[1] == 42.2081 );
// }
// {
//  namespace lapack = multi::lapack;
//  multi::array<double, 2> A = {
//      {167.413, 126.804, 125.114},
//      {NAN    , 167.381, 126.746},
//      {NAN    , NAN    , 167.231}
//  };
//  auto W = lapack::syev(multi::blas::filling::upper, A);
//  BOOST_TEST( A[2][1] == -0.579092 );
//  BOOST_TEST( W[1] == 42.2081 );
// }
// {
//  multi::array<double, 2> const A = {
//      {167.413, 126.804, 125.114},
//      {NAN    , 167.381, 126.746},
//      {NAN    , NAN    , 167.231}
//  };
//  multi::array<double, 1> W(size(A));
//  namespace lapack = multi::lapack;
//  auto A_copy = lapack::syev(lapack::filling::upper, A, W);
//  BOOST_TEST( A[1][2] == 126.746 );
//  BOOST_TEST( A_copy[2][1] == -0.579092 );
//  BOOST_TEST( W[1] == 42.2081 );
// }
// {
//  multi::array<double, 2> A = {
//      {167.413, 126.804, 0.},
//      {NAN    , 167.381, 0.},
//      {NAN    , NAN    , 0.}
//  };
//  multi::array<double, 1> W(size(A));
//  namespace lapack = multi::lapack;
//  auto&& A_ref = lapack::syev(lapack::filling::upper, A, W);
//  BOOST_TEST( size(A_ref)==3 );
//  BOOST_TEST( W[0]==0. );
// }
// {
//  multi::array<double, 2> A = {
//      {1. , 1.,  1.},
//      {NAN, 2 ,  1.},
//      {NAN, NAN, 1.}
//  };
//  multi::array<double, 1> W(size(A));
//  namespace lapack = multi::lapack;
//  auto&& A_ref = lapack::syev(lapack::filling::upper, A, W);
//  print(A_ref);
//  BOOST_TEST( size(A_ref)==3 );
//  BOOST_TEST( W[0]==0. );
// }
// {
//  multi::array<double, 2> A = {{5.}};
//  multi::array<double, 1> W(size(A));
//  namespace lapack = multi::lapack;
//  lapack::syev(lapack::filling::upper, A, W);
//  BOOST_TEST( A[0][0] == 1. );
//  BOOST_TEST( W[0]==5. );
// }
// {
//  namespace lapack = multi::lapack;
//  multi::array<double, 2> A;
//  multi::array<double, 1> W(size(A));
//  lapack::syev(lapack::filling::upper, A, W);
// }
// {
//  multi::array<double, 2> const A = {
//      {167.413, 126.804, 125.114},
//      {NAN    , 167.381, 126.746},
//      {NAN    , NAN    , 167.231}
//  };
//  multi::array<double, 1> W(size(A));
//  namespace lapack = multi::lapack;
//  auto sys = lapack::syev(lapack::filling::upper, A);
//  BOOST_TEST( A[1][2] == 126.746 );
//  BOOST_TEST( sys.eigenvectors[2][1] == -0.579092 );
//  BOOST_TEST( sys.eigenvalues[1] == 42.2081 );
// }
// {
//  multi::array<double, 2> const A = {
//      {167.413, 126.804, 125.114},
//      {NAN    , 167.381, 126.746},
//      {NAN    , NAN    , 167.231}
//  };
//  multi::array<double, 1> W(size(A));
//  namespace lapack = multi::lapack;
//  auto [eigenvecs, eigenvals] = lapack::syev(lapack::filling::upper, A);
//  BOOST_TEST( A[1][2] == 126.746 );
//  BOOST_TEST( eigenvecs[2][1] == -0.579092 );
//  BOOST_TEST( eigenvals[1] == 42.2081 );
// }

// #endif
