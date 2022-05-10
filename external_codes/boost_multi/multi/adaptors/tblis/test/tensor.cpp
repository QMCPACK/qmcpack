// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;-*-
// © Alfredo A. Correa 2019-2021

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi BLAS dot"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include "../../tblis.hpp"

namespace multi = boost::multi;

BOOST_AUTO_TEST_CASE(blis_matrix){

	namespace tblis = multi::tblis;
	using namespace multi::tblis;

	auto const A = []{
		multi::array<double, 2> _({5, 2}); std::iota(_.elements().begin(), _.elements().end(), 0.);
		return _;
	}();

	auto const B = []{
		multi::array<double, 2> _({2, 7}); std::iota(_.elements().begin(), _.elements().end(), 0.);
		return _;
	}();

	// now the check
	multi::array<double, 2> C_gold({5, 7}, 0.);
	
	assert( extension(C_gold) == extension(A) );
	assert( extension(C_gold[0]) == extension(B[0]) );
	assert( extension(B) == extension(A[0]) );
	for(auto a : extension(C_gold)){
		for(auto b : extension(C_gold[0])){
			for(auto c : extension(B)){
				C_gold[a][b] += A[a][c]*B[c][b];
			}
		}
	}

	{
		multi::array<double, 2> C({5, 7}, 0.);
		// C[abcd] += A[cebf] B[afed]
		tblis::mult(tblis::matrix(A), tblis::matrix(B), tblis::matrix(C));
		BOOST_REQUIRE( C_gold == C );
	}
	{
		multi::array<double, 2> C({5, 7}, 0.);
		tblis::mult(tblis::tensor(A), "ac", tblis::tensor(B), "cb", tblis::tensor(C), "ab");
		BOOST_REQUIRE( C_gold == C );
	}
	{
		multi::array<double, 2> C({5, 7}, 0.);
		tblis::mult(tblis::tensor(A)["ac"], tblis::tensor(B)["cb"], tblis::tensor(C)["ab"]);
		BOOST_REQUIRE( C_gold == C );
	}
	{
		multi::array<double, 2> C({5, 7}, 0.);
		tblis::mult(tblis::tensor(A)["ac"], tblis::tensor(B)["cb"], tblis::tensor(C)["ab"]);
		BOOST_REQUIRE( C_gold == C );
	}
	{
		multi::array<double, 2> C({5, 7}, 0.);
		using namespace tblis::indices;
		tblis::mult(tblis::tensor(A)(a, c), tblis::tensor(B)(c, b), tblis::tensor(C)(a, b));
		BOOST_REQUIRE( C_gold == C );
	}
	{
		multi::array<double, 2> C({5, 7}, 0.);
		using namespace tblis::indices;
		tblis::mult(A(a, c), B(c, b), C(a, b));
	//	BOOST_REQUIRE( C_gold == C );
	}
}

BOOST_AUTO_TEST_CASE(tblis_tensor){
	namespace multi = boost::multi;
	namespace tblis = multi::tblis;

	auto const A = []{
		multi::array<double, 4> A({2, 5, 10, 9});
		std::iota(A.data_elements(), A.data_elements() + A.num_elements(), 0.);
		return A;
	}();

	auto const B = []{
		multi::array<double, 4> B({8, 9, 5, 7});
		std::iota(B.data_elements(), B.data_elements() + B.num_elements(), 0.);
		return B;
	}();

	auto const C_gold = [&A, &B]{
		multi::array<double, 4> _({8, 10, 2, 7}, 0.);
		// this computers C_check[abcd] += A[cebf] B[afed]
		for(auto a = 0; a != 8; ++a){
			for(auto b = 0; b != 10; ++b){
				for(auto c = 0; c != 2; ++c){
					for(auto d = 0; d != 7; ++d){

						for(auto e = 0; e != 5; ++e){
							for(auto f = 0; f != 9; ++f){
								_[a][b][c][d] += A[c][e][b][f]*B[a][f][e][d];
							}
						}

					}
				}
			}
		}
		return _;
	}();

	{
		multi::array<double, 4> C({8, 10, 2, 7}, 0.);
		{
			using namespace tblis::indices;
			tblis::mult( A(c, e, b, f), B(a, f, e, d), C(a, b, c, d) );
		}
		BOOST_REQUIRE( C_gold == C );
	}
	#if defined(__clang__)
	{
		multi::array<double, 4> C({8, 10, 2, 7}, 0.);
		{
			using namespace tblis::indices::greek;
			tblis::mult( A(γ, ε, β, ζ), B(α, ζ, ε, δ), C(α, β, γ, δ) );
		}
		BOOST_REQUIRE( C_gold == C );
	}
	#endif
	{
		multi::array<double, 4> C({8, 10, 2, 7}, 0.);
		{
			using namespace tblis::indices;
			tblis::mult( tblis::tensor(A)(c, e, b, f), tblis::tensor(B)(a, f, e, d), tblis::tensor(C)(a, b, c, d) );
		}
		BOOST_REQUIRE( C_gold == C );
	}


}

