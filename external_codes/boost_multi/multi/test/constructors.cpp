#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4-*-
echo $X
$CXXX $CXXFLAGS $0 -o $0.$X -lboost_unit_test_framework&&$0.$X&&rm $0.$X;exit
#endif
// © Alfredo A. Correa 2019-2020

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi constructors"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include "../array.hpp"

#include<complex>
#include<functional>
#include<iostream>
#include<numeric>
#include<vector>

namespace multi = boost::multi;

using complex = std::complex<double>;

struct multiplies_bind1st{
	explicit multiplies_bind1st(multi::array<complex, 2>&& m) : m_(std::move(m)){} // this produces a bug in nvcc11.0
private:
	multi::array<complex, 2> m_;
};

BOOST_AUTO_TEST_CASE(multi_construct_1d){
	multi::static_array<double, 1> A(multi::array<double, 1>::extensions_type{10}, 1.);
	BOOST_REQUIRE( size(A) == 10 );
	BOOST_REQUIRE( A[1] == 1. );
}

BOOST_AUTO_TEST_CASE(multi_constructors_inqnvcc_bug){
	multi::array<complex, 2> m({10, 10});
	multiplies_bind1st(std::move(m));
}

BOOST_AUTO_TEST_CASE(multi_constructors_1d){
	{
		multi::array<double, 1> A(10); 
		BOOST_REQUIRE( size(A)==10 );
	}
	{
		multi::array<double, 1> A(10, double{}); 
		BOOST_REQUIRE( size(A)==10 );
		BOOST_REQUIRE( A[5]== double{} );
	}
	{
		multi::array<double, 1> A(10, double{}); 
		BOOST_REQUIRE( size(A)==10 );
		BOOST_REQUIRE( A[5]== double{} );
	}
	{
	#if defined(__cpp_deduction_guides)
		multi::array A(10, double{}); 
		BOOST_REQUIRE( size(A)==10 );
		BOOST_REQUIRE( A[5]== double{} );
	#endif
	}
}

BOOST_AUTO_TEST_CASE(multi_constructors){
{//multi::array<double, 1> A({10}); assert(size(A)==1); // warning in clang
}{//multi::array<double, 1> A({10}, double{}); assert(size(A)==10); // warning in clang
}{//multi::array<double, 1> A({10}, double{}); assert(size(A)==10); // warning in clang
}{//multi::array<double, 1> A({10}, 0.); assert(size(A)==10); // warning in clang
}{//multi::array<double, 1> A({10}, {}); assert(size(A)==10); // error ambiguous 
}{ multi::array<std::size_t, 1> A = {10}    ; BOOST_REQUIRE( size(A)==1 and A[0]==10 );
}{ multi::array<int        , 1> A = {10}    ; BOOST_REQUIRE( size(A)==1 and A[0]==10 );
}{ multi::array<double     , 1> A = {10}    ; BOOST_REQUIRE( size(A)==1 and A[0]==10 );
}{ multi::array<std::size_t, 1> A({10})     ; BOOST_REQUIRE( size(A)==1 and A[0]==10 );
}{ multi::array<int        , 1> A({10})     ; BOOST_REQUIRE( size(A)==1 and A[0]==10 );
}{ multi::array<double     , 1> A({10})     ; BOOST_REQUIRE( size(A)==1 and A[0]==10 );
//}{ multi::array<std::size_t, 1> A({{10}})   ; assert( size(A)==1 and A[0]==10 );  // clang warns about double bracked
//}{ multi::array<int        , 1> A({{10}})   ; assert( size(A)==1 and A[0]==10 );  // clang warns about double bracked
//}{ multi::array<double     , 1> A({{10}})   ; assert( size(A)==1 and A[0]==10 );  // clang warns about double bracked
}{ multi::array<std::size_t, 1> A({0, 10})  ; BOOST_REQUIRE( size(A)==2 );
}{ multi::array<int        , 1> A({0, 10})  ; BOOST_REQUIRE( size(A)==2 );
}{ multi::array<double     , 1> A({0, 10})  ; BOOST_REQUIRE( size(A)==2 );
}
}

//}{ multi::array<std::size_t, 1> A({{0, 10}}); assert( size(A)==2 ); // ambiguous in gcc (error)
//}{ multi::array<int        , 1> A({{0, 10}}); assert( size(A)==2 ); // ambiguous in gcc (error)
//}{ multi::array<double     , 1> A({{0, 10}}); assert( size(A)==2 ); // ambiguous in gcc (error)

//}{ multi::array<std::size_t, 2> A = {{1, 2}, {3, 4}}; assert( size(A)==2 and A[0][0]==1 );
//}{ multi::array<int,         2> A = {{1, 2}, {3, 4}}; assert( size(A)==2 and A[0][0]==1 );
//}{ multi::array<double,      2> A = {{1, 2}, {3, 4}}; assert( size(A)==2 and A[0][0]==1 );
//}{ multi::array<std::size_t, 2> A({{1, 2}, {3, 4}}) ; assert( size(A)==2 and A[0][0]==1 );
//}{ multi::array<int,         2> A({{1, 2}, {3, 4}}) ; assert( size(A)==2 and A[0][0]==1 );
//}{ multi::array<double,      2> A({{1, 2}, {3, 4}}) ; assert( size(A)==2 and A[0][0]==1 );
//}
//{
//	multi::array<double, 3> A(multi::index_extensions<3>{4, 5, 6}, 8.);
//	multi::array<double, 3> B(multi::iextensions<3>{4, 5, 6}, 8.);
//	multi::array<double, 3> B2(multi::iextensions<3>{std::array<multi::size_type, 3>{4, 5, 6}}, 8.);
//	multi::array<double, 3> B3(multi::iextensions<3>{std::array<multi::size_type, 3>{4, 5, 6}});
//	multi::array<double, 3> C({4, 5, 6}, 8.);
//	multi::array<double, 3> D(std::array<multi::size_type, 3>{4, 5, 6}, 8.);
//	multi::array<double, 3> D1(std::tuple<multi::size_type, multi::size_type, multi::size_type>{4, 5, 6}, 8.);
//	multi::array<double, 3> D2(std::tuple<multi::size_type, multi::size_type, multi::size_type>{4, 5, 6});

//	multi::array<double, 3> E({4, 5, 6}); assert(size(E)==4);
//	multi::array<double, 3> E2({4, 5, 6}); assert(size(E)==4);

//}
//{
//	multi::array<double, 2> A(multi::index_extensions<2>{8, 8}, 8.);
//	assert( size(A) == 8 );
//	assert( std::get<0>(sizes(A)) == 8 );
//	assert( std::get<1>(sizes(A)) == 8 );

//}{
//	multi::array<double, 2> A({8, 8}, 8.);
//	assert( size(A) == 8 );
//	assert( std::get<0>(sizes(A)) == 8 );
//	assert( std::get<1>(sizes(A)) == 8 );
//}
// {  multi::static_array<double, 1> A     ; assert( is_empty(A) );
//}{  multi::static_array<double, 1> A{}   ; assert( is_empty(A) );
//}{  multi::static_array<double, 1> A = {}; assert( is_empty(A) ); 
//}{  multi::static_array<double, 2> A     ; assert( is_empty(A) );
//}{  multi::static_array<double, 2> A{}   ; assert( is_empty(A) );
//}{  multi::static_array<double, 2> A = {}; assert( is_empty(A) );
//}{  multi::static_array<double, 3> A     ; assert( is_empty(A) );
//}{  multi::static_array<double, 3> A{}   ; assert( is_empty(A) );
//}{  multi::static_array<double, 3> A = {}; assert( is_empty(A) );
//}{  multi::static_array<double, 3> A, B  ; assert( A == B );
//}{
//	multi::array<double, 1> A1 = {0.0, 1.0, };
//	assert( size(A1) == 2 );
//	assert( A1[1] == 1.0 );
//}

//{
//	double a[4][5];
//	double const b[4][5] = {
//		{ 0,  1,  2,  3,  4}, 
//		{ 5,  6,  7,  8,  9}, 
//		{10, 11, 12, 13, 14}, 
//		{15, 16, 17, 18, 19}
//	};;
//	multi::array_ptr<double, 2> ap(&a[0][0], {4, 5});
//	multi::array_ptr
//#if not defined(__cpp_deduction_guides)
//		<double, 2, double const*>
//#endif
//		bp(&b[0][0], {4, 5})
//	;
//	*ap = *bp;
//	BOOST_REQUIRE( (*ap).rotated() == (*bp).rotated() );
//}
//{
//	multi::array<double, 2> A = {{0, 3}, {0, 5}};
//	assert( size(A) == 2 );
//	assert( size(A[0]) == 2 );
//	assert( A[1][1] == 5 );
//}
//{
//	multi::array<double, 2> A({{0, 3}, {0, 5}});
//	assert( size(A) == 2 );
//	assert( size(A[0]) == 2 );
//}
////{
////	std::vector<multi::static_array<double, 2>> v(9, {multi::index_extensions<2>{8, 8}});
////	#if __cpp_deduction_guides
////	std::vector w(9, multi::static_array<double, 2>({8, 8}));
////	#endif
////}
// {  multi::array<double, 1, std::allocator<double>> A{std::allocator<double>{}}; assert( is_empty(A) );
//}{  multi::array<double, 2, std::allocator<double>> A{std::allocator<double>{}}; assert( is_empty(A) );
//}{  multi::array<double, 3, std::allocator<double>> A{std::allocator<double>{}}; assert( is_empty(A) );
//}{ multi::array<double, 1> A(3, {});  assert( size(A)==3 );
//}{ multi::array<double, 1> A(3, 99.); assert( size(A)==3 and A[2]==99. );
//}{ multi::array<double, 1> A({3});    assert( size(A)==1 and A[0]==3 );
//}{// multi::array<double, 1> A({{3}});  assert( size(A)==1 and A[0]==3 );
//}{ multi::array<double, 1> A(multi::iextensions<1>{3}); assert(size(A)==3);
//}{ multi::array<double, 1> A(multi::array<double, 1>::extensions_type{3}); assert(size(A)==3);
//}

#if 0


//#if not defined(__INTEL_COMPILER)
}{	multi::array<double, 1> A({3}); assert( size(A)==1 and A[0]==3. );  // uses init_list
}{	multi::array<double, 1> A({3}); assert( size(A)==1 and A[0]==3. );  // uses init_list
//#endif
//#if not defined(__INTEL_COMPILER)
}{  multi::array<double, 1> A({3.}); assert( size(A)==1 and A[0]==3. ); // uses init_list
}{  multi::array<double, 1> A = {3l}; assert( size(A)==1 and A[0]==3. ); // uses init_list
//#else
}{  multi::array<double, 1> A(3l, 0.); assert( size(A)==3 and A[0]==0. ); // gives warning in clang++ A({3l}, 0.);
//#endif
}{  multi::array<double, 1> A(multi::index_extensions<1>{{0, 3}}); assert( size(A)==3 and A[0]==0 );
#if (!defined(__INTEL_COMPILER)) && (defined(__GNUC) && __GNU_VERSION__ >= 600)
//}{  multi::array<double, 1> A({{0l, 3l}}); cout<<size(A)<<std::endl; assert( size(A)==3 and A[1]==0. ); //uses init_list
#endif
}{  // multi::array<double, 1, std::allocator<double>> A(multi::index_extensions<1>{2}, std::allocator<double>{}); assert( size(A)==2 );
}{  // multi::array<double, 1, std::allocator<double>> A(multi::index_extensions<1>{{0, 3}}, std::allocator<double>{}); assert( size(A)==3 );
}{  // multi::array<double, 1, std::allocator<double>> A(multi::iextensions<1>{2}, std::allocator<double>{}); assert( size(A)==2 );
}{  // multi::array<double, 1, std::allocator<double>> A(multi::iextensions<1>{{0, 3}}, std::allocator<double>{}); assert( size(A)==3 );
#if not defined(__INTEL_COMPILER) and (defined(__GNUC) and __GNU_VERSION >= 600)
}{  multi::array<double, 2> A({2, 3}); assert( num_elements(A)==6 );
#endif
}{  multi::array<double, 2> A(multi::iextensions<2>{2, 3}); assert( num_elements(A)==6 );
//}{  multi::array<double, 2> A({2, 3}); assert( num_elements(A)==6 and size(A)==2 and std::get<1>(sizes(A))==3 );
}{  multi::array<double, 2> A(multi::index_extensions<2>{{0,2}, {0,3}}); assert( num_elements(A)==6 );
#if not defined(__INTEL_COMPILER) and (defined(__GNUC__) and __GNU_VERSION__ >= 600)
}{  multi::array<double, 2, std::allocator<double>> A({2, 3}, std::allocator<double>{}); assert( num_elements(A)==6 );
#endif
}{  multi::array<double, 2, std::allocator<double>> A(multi::iextensions{2, 3}, std::allocator<double>{}); assert( num_elements(A)==6 );
#if not defined(__INTEL_COMPILER) and (defined(__GNUC__) and (__GNU_VERSION >= 600))
}{  multi::array<double, 3> A({2, 3, 4}); assert( num_elements(A)==24 and A[1][2][3]==0 );
#endif
}{  multi::array<double, 3> A(multi::iextensions<3>{2, 3, 4}); assert( num_elements(A)==24 and A[1][2][3]==0 );
#if not defined(__INTEL_COMPILER) and (defined(__GNUC__) and __GNU_VERSION__ >= 600 )
}{  multi::array<double, 3> A({{0, 2}, {0, 3}, {0, 4}}); assert( num_elements(A)==24 and A[1][2][3]==0 );
#endif
}{  multi::array<double, 3> A(multi::iextensions<3>{{0, 2}, {0, 3}, {0, 4}}); assert( num_elements(A)==24 and A[1][2][3]==0 );
#if (not defined(__INTEL_COMPILER)) and (defined(__GNUC__) and __GNU_VERSION__ >= 600)
}{  multi::array<double, 3, std::allocator<double>> A({2, 3, 4}, std::allocator<double>{}); assert( num_elements(A)==24 );
#endif
}

 {  multi::array<double, 1> A(multi::iextensions<1>{3}, 3.1); assert( size(A)==3 and A[1]==3.1 );
}{//multi::array<double, 1> A({3}, 3.1); assert( size(A)==3 and A[1]==3.1 ); // warning in clang
}{//multi::array<double, 1> A({3l}, 3.1); assert( size(A)==3 and A[1]==3.1 ); // warning in clang
}{  multi::array<double, 1> A( {{0,3}}, 3.1); assert( size(A)==3 and A[1]==3.1 );
}{  multi::array<double, 1> A(multi::iextension(3), 3.1); assert( size(A)==3 and A[1]==3.1 );
}{  multi::array<double, 1> A(3l, 3.1); assert( size(A)==3 and A[1]==3.1 );
}{  multi::array<double, 1> A(3, 3.1); assert( size(A)==3 and A[1]==3.1 );
}{  multi::array<double, 1> A({{0, 3}}, 3.1); assert( size(A)==3 and A[1]==3.1 );
#if (not defined(__INTEL_COMPILER)) and (defined(__GNUC__) and __GNU_VERSION__ >=600)
}{  multi::array<double, 2> A({2, 3}, 3.1); assert( num_elements(A)==6 and A[1][2]==3.1 );
#endif
}{  multi::array<double, 2> A(multi::iextensions<2>{2, 3}, 3.1); assert( num_elements(A)==6 and A[1][2]==3.1 );
#if (not defined(__INTEL_COMPILER)) and (defined(__GNUC__) and __GNU_VERSION__ >=600)
}{  multi::array<double, 2> A({{0,2}, {0,3}}, 3.1); assert( num_elements(A)==6 and A[1][2]==3.1 );
#endif
}{  multi::array<double, 2> A(multi::iextensions<2>{{0,2}, {0,3}}, 3.1); assert( num_elements(A)==6 and A[1][2]==3.1 );
#if (not defined(__INTEL_COMPILER)) and (defined(__GNUC__) and __GNU_VERSION__ >=600)
}{  multi::array<double, 3> A({2, 3, 4}, 3.1); assert( num_elements(A)==24 and A[1][2][3]==3.1 );
#endif
}{  multi::array<double, 3> A(multi::iextensions<3>{2, 3, 4}, 3.1); assert( num_elements(A)==24 and A[1][2][3]==3.1 );
}
#endif
//}

