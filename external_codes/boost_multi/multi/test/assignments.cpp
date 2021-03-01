#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;-*-
$CXX $CXXFLAGS $0 -o $0.$X -lboost_unit_test_framework&&$0.$X&&rm $0.$X;exit
#endif
// © Alfredo A. Correa 2019-2020

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi assignments"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include<complex>

#include "../array.hpp"

namespace multi = boost::multi;

multi::array_ref<double, 2> make_ref(double* p){return {p, {5, 7}};}

BOOST_AUTO_TEST_CASE(equality_1D){
	multi::array<double, 1> A = {1., 2., 3.};
	multi::array<double, 1> B = {1., 2., 3.};
	BOOST_REQUIRE( A == B );
	BOOST_REQUIRE( not (A != B) );

	BOOST_REQUIRE( A() == B() );
	BOOST_REQUIRE( not (A() != B()) );
}

BOOST_AUTO_TEST_CASE(equality_2D){
	multi::array<double, 2> A = {
		{1., 2., 3.},
		{4., 5., 6.}
	};
	multi::array<double, 2> B = {
		{1., 2., 3.},
		{4., 5., 6.}
	};
	BOOST_REQUIRE( A == B );
	BOOST_REQUIRE( not (A != B) );

	BOOST_REQUIRE( A() == B() );
	BOOST_REQUIRE( not (A() != B()) );
	
	BOOST_REQUIRE( A[0] == B[0] );
	BOOST_REQUIRE( not (A[0] != B[0]) );
}

BOOST_AUTO_TEST_CASE(multi_copy_move){
	multi::array<double, 2> A({3, 3}, 0.);
	multi::array<double, 2> B = A;//identy();//multi::array<double, 2>({3, 3}, 0.);// = eye<double>({3, 
	BOOST_REQUIRE( A == B );

	auto A_data = A.data_elements();
	multi::array<double, 2> C = std::move(A);
	BOOST_REQUIRE( is_empty(A) );
	BOOST_REQUIRE( A_data = C.data_elements() );
	
	multi::array<double, 2> D(std::move(B));
	BOOST_REQUIRE( is_empty(B) );
}

#if 1
BOOST_AUTO_TEST_CASE(range_assignment){
{
	auto r = multi::make_extension_t(10l);
	multi::array<double, 1> v(r.begin(), r.end());
	BOOST_REQUIRE( r.size() == v.size() );
	BOOST_REQUIRE( v[1] = 10 );
}
{
	multi::array<double, 1> v(10);
	auto r = extension(v);
	v.assign(r.begin(), r.end());
	BOOST_REQUIRE( v[1] == 1 );
}
}

BOOST_AUTO_TEST_CASE(rearranged_assignment){
	multi::array<double, 4> tmp({14, 14, 7, 4});
	multi::array<double, 5> src({2, 14, 14, 7, 2}); src[0][1][2][3][1] = 99.;

	BOOST_REQUIRE( extensions(tmp.unrotated().partitioned(2).transposed().rotated()) == extensions(src) );
}

BOOST_AUTO_TEST_CASE(rvalue_assignments){
	using complex = std::complex<double>;

	std::vector<double> const v1(200, 99.);
	std::vector<complex> v2(200);
	auto linear1 = [&]{return multi::array_cptr<double, 1>(v1.data(), 200);};
	auto linear2 = [&]{return multi::array_ptr<complex, 1>(v2.data(), 200);};
	*linear2() = *linear1();

}

#if 0 // self-move-assigment is a standard warning in clang (-Wmove)
BOOST_AUTO_TEST_CASE(self_assigment){
	multi::array<double, 1> A = {1., 2., 3.};
	A = std::move(A);
	std::cout << A[0] << std::endl;
	BOOST_REQUIRE( A.empty() );

	multi::array<double, 2> B = {{1., 2., 3.},{2.,3.,4.}};
	B = std::move(B);
	BOOST_REQUIRE( B.empty() );
}
#endif

BOOST_AUTO_TEST_CASE(assignments){
	{
		std::vector<double> v(5*7, 99.);
		constexpr double val = 33.;
		multi::array<double, 2> A({5, 7}, val);
		multi::array_ref<double, 2>(v.data(), {5, 7}) = A;
		BOOST_REQUIRE( v[9] == val );
		BOOST_REQUIRE( not v.empty() );
		BOOST_REQUIRE( not is_empty(A) );

		multi::array<double, 1> V;
		BOOST_REQUIRE( V.empty() );
	}
	{
		std::vector<double> v(5*7, 99.), w(5*7, 33.);

		multi::array_ptr<double, 2> Bp{w.data(), {5, 7}};
		make_ref(v.data()) = *Bp;
		make_ref(v.data()) = Bp->sliced(0, 5);

		BOOST_REQUIRE( v[9] == 33. );
	}
	{
		std::vector<double> v(5*7, 99.), w(5*7, 33.);

		make_ref(v.data()) = make_ref(w.data());

		BOOST_REQUIRE( v[9] == 33. );
	}
}

template<class T, class Allocator = std::allocator<T> > // T must be specified, double, complex<double>
multi::array<T, 2, Allocator> eye(multi::iextensions<2> ie, Allocator alloc = {}){
	multi::array<T, 2, Allocator> ret(ie, 0., alloc);
	ret.diagonal().fill(1.);
	return ret;
}

BOOST_AUTO_TEST_CASE(assigment_temporary){
	multi::array<double, 2> Id = eye<double>({3, 3});
	BOOST_REQUIRE( Id == eye<double>({3, 3}) );
	BOOST_REQUIRE( Id[1][1] == 1 );
	BOOST_REQUIRE( Id[1][0] == 0 );
}

#endif

