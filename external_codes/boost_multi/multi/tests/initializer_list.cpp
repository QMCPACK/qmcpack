#ifdef COMPILATION_INSTRUCTIONS
$CXX -O3 -std=c++17 -Wall -Wextra -Wpedantic $0 -o $0.x && $0.x $@ && rm -f $0.x; exit
#endif

#include "../array.hpp"

#include<complex>
#include<iostream>

namespace multi = boost::multi;

int main(){
#if 1
{
	multi::array<double, 1> const A = {1.2,3.4,5.6}; 
	assert( size(A) == 3 );
	assert( A[2] == 5.6 );
}
{
	double const a[3] = {1.1, 2.2, 3.3};
	using multi::num_elements;
	assert( num_elements(a) == 3 );
	multi::array<double, 1> const A(a);
	assert(size(A) == 3);
}
{
//	multi::array<double, 1> const A((double[]){1.1,2.2,3.3}); // warning: ISO C++ forbids compound-literals [-Wpedantic]
//	assert( size(A)==3 and A[1] == 2.2 );
}
{
//	std::array<double, 3> a = {1.1,2.2,3.3};
//	multi::array<double, 1> const A(a); assert( size(A)==3 and A[1]==2.2 );
}
{
//	multi::array<double, 1> const A(std::array<double, 3>{1.1,2.2,3.3}); assert( size(A)==3 and A[1]==2.2 );
}
{
	multi::array<double, 2> const A = {
		{ 1.2,  2.4, 3.6},
		{11.2, 34.4, 5.6},
		{15.2, 32.4, 5.6}
	};
	multi::array<double, 2> B; B = A;
	assert( size(A) == 3 );
	assert( size(A) == 3 and size(A[0]) == 3 );
	assert( A[1][1] == 34.4 );
}
{
	multi::array<double, 2> const A = {
		{ 1.2,  2.4},
		{11.2, 34.4},
		{15.2, 32.4}
	};
	assert( size(A) == 3 and size(A[0]) == 2 );
	assert( A[1][1] == 34.4 );
}
{
	double const a[3][2] = {
		{ 1.2,  2.4},
		{11.2, 34.4},
		{15.2, 32.4}
	};
	multi::array<double, 2> A(a);
}
{
	double const staticA[][2] = 
	{	{ 1.2,  2.4},
		{11.2, 34.4},
		{15.2, 32.4}};
	multi::array<double, 2> A = staticA;
	assert( size(A) == 3 );
}
{
	multi::array<double, 2> A = 
		(double const[][2]) // warns with -Wpedantic
	{	{ 1.2,  2.4},
		{11.2, 34.4},
		{15.2, 32.4}};
	assert( size(A) == 3 );
}
{
	multi::array<double, 2> A = 
		#if defined(__INTEL_COMPILER)
		(double const[3][4])
		#endif
	{	{ 1.2,  2.4},
		{11.2, 34.4},
		{15.2, 32.4}};
}
{
//	std::array<std::array<double, 2>, 3> a = {{
//		{{1.,2.}},
//		{{2.,4.}},
//		{{3.,6.}}
//	}};
//	multi::array<double, 2> A(a);
//	assert( num_elements(A) == 6 and A[2][1] == 6. );
}
{
//	multi::array<double, 2> A(
//		std::array<std::array<double, 2>, 3>{{
//			{{1.,2.}},
//			{{2.,4.}},
//			{{3.,6.}}
//		}}
//	);
}
{
	multi::array<double, 3> const A = {  // warning: ISO C++ forbids compound-literals [-Wpedantic]
		{{ 1.2,  0.}, { 2.4, 1.}},
		{{11.2,  3.}, {34.4, 4.}},
		{{15.2, 99.}, {32.4, 2.}}
	};
	assert( A[1][1][0] == 34.4 and A[1][1][1] == 4.   );
}
#if __cpp_deduction_guides
{
	multi::array const A = {1.2,3.4,5.6}; 
	assert( size(A) == 3 );
	assert( A[2] == 5.6 );
}
{
	multi::array const A = {
		{ 1.2,  2.4, 3.6},
		{11.2, 34.4, 5.6},
		{15.2, 32.4, 5.6}
	};
	assert( size(A) == 3 );
	assert( size(A) == 3 and size(A[0]) == 3 );
	assert( A[1][1] == 34.4 );
}
{
	multi::array const A = {
		{ 1.2,  2.4},
		{11.2, 34.4},
		{15.2, 32.4}
	};
	assert( size(A) == 3 and size(A[0]) == 2 );
	assert( A[1][1] == 34.4 );
}
{
	multi::array const A = {
		{{ 1.2,  0.}, { 2.4, 1.}},
		{{11.2,  3.}, {34.4, 4.}},
		{{15.2, 99.}, {32.4, 2.}}
	};
	assert( A[1][1][0] == 34.4 and A[1][1][1] == 4.   );
}
#endif

{
	std::complex<double> const I(0.,1.);
	multi::array<std::complex<double>, 2> b = {
		{2. + 1.*I, 1. + 3.*I, 1. + 7.*I},
		{3. + 4.*I, 4. + 2.*I, 0. + 0.*I}
	};
}
{	
	multi::array<double, 1> A1 = {1., 2., 3.}; 
	assert(num_elements(A1)==3 and A1[1]==2.);
	multi::array<double, 1> B1 = {3., 4.}; assert(num_elements(B1)==2 and B1[1]==4.);
	multi::array<double, 1> C1 = {0, 4}; assert(num_elements(C1)==2 and C1[1]==4.);
	multi::array<double, 1> D1 = {0l, 4l}; assert(num_elements(D1)==2 and D1[1]==4.);
//	multi::array<double, 1> E1 = {{0, 4}}; assert(num_elements(E1)==2 and E1[1]==4.); // [X] icc19
	multi::array<double, 1> F1({0, 4}); assert(num_elements(F1)==2 and F1[1]==4.);
	multi::array<double, 2> A2 = {
		{1., 2., 3.}, 
		{4., 5., 6.}
	}; 
	assert(num_elements(A2)==6 and A2[1][1]==5.);
	multi::array<double, 3> A3 = {
		{ { 1.,  2.,  3.}, 
		  { 4.,  5.,  6.} },
		{ { 7.,  8.,  9.}, 
		  {10., 11., 12.} }
	};
	assert(num_elements(A3)==12 and A3[1][1][1]==11.);
	
	multi::array<std::string, 3> B3 = {
		{ {"000", "001", "002"}, 
		  {"010", "011", "012"} },
		{ {"100", "101", "102"}, 
		  {"110", "111", "112"} }
	};
	assert( num_elements(B3)==12 and B3[1][0][1] == "101" );
}
#if __cpp_deduction_guides
 {	multi::array A = {1., 2., 3.}; assert( multi::rank<decltype(A)>{}==1 and num_elements(A)==3 and A[1]==2. ); static_assert( typename decltype(A)::rank{}==1 );
}{	multi::array A = {1., 2.};     assert( multi::rank<decltype(A)>{}==1 and num_elements(A)==2 and A[1]==2. ); assert( multi::rank<decltype(A)>{}==1 );
}{	multi::array A = {0, 2};       assert( multi::rank<decltype(A)>{}==1 and num_elements(A)==2 and A[1]==2. ); assert( multi::rank<decltype(A)>{}==1 );
}{	multi::array A = {9.};         assert( multi::rank<decltype(A)>{}==1 and num_elements(A)==1 and A[0]==9. ); assert( multi::rank<decltype(A)>{}==1 );
}{	multi::array A = {9};          assert( multi::rank<decltype(A)>{}==1 and num_elements(A)==1 and A[0]==9. ); assert( multi::rank<decltype(A)>{}==1 );
}{	multi::array A = {
		{1., 2., 3.}, 
		{4., 5., 6.}
	}; assert( multi::rank<decltype(A)>{}==2 and num_elements(A)==6 );
}{
	multi::array A = {1., 2., 3.};
	multi::array B = {1., 2.};
//	multi::array C = {A, A, B, A}; assert( dimensionality(C) == 1 and num_elements(C) == 3 );
}{
	multi::array B3 = {
		{ {"000", "001", "002"}, 
		  {"010", "011", "012"} },
		{ {"100", "101", "102"}, 
		  {"110", "111", "112"} }
	};
	static_assert( std::is_same<decltype(B3)::element, char const*>{}, "!");
	static_assert( not std::is_same<decltype(B3)::element, std::string>{}, "!");
	auto C3 = B3; 
}
#endif
#endif

}

