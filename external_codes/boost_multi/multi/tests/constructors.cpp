#ifdef COMPILATION_INSTRUCTIONS
c++ -O3 -std=c++17 -Wall -Wextra -Wpedantic `#-Wfatal-errors` $0 -o $0.x && $0.x $@ &&rm -f $0.x; exit
#endif

#include<iostream>
#include "../array.hpp"

#include<queue>
#include<vector>
#include<numeric> // iota

namespace multi = boost::multi;
using std::cout;

int main(){

using boost::multi::size;

 {	multi::array<double, 1> A     ; assert( empty(A) );
}{	multi::array<double, 1> A{}   ; assert( empty(A) );
}{	multi::array<double, 1> A = {}; assert( empty(A) );
}{  multi::array<double, 2> A     ; assert( empty(A) );
}{  multi::array<double, 2> A{}   ; assert( empty(A) );
}{  multi::array<double, 2> A = {}; assert( empty(A) );
}{  multi::array<double, 3> A     ; assert( empty(A) );
}{  multi::array<double, 3> A{}   ; assert( empty(A) );
}{  multi::array<double, 3> A = {}; assert( empty(A) );
}{  multi::array<double, 3> A, B  ; assert( A == B );
}

 {	multi::array<double, 1, std::allocator<double>> A{std::allocator<double>{}}; assert( empty(A) );
}{	multi::array<double, 2, std::allocator<double>> A{std::allocator<double>{}}; assert( empty(A) );
}{	multi::array<double, 3, std::allocator<double>> A{std::allocator<double>{}}; assert( empty(A) );
}

 {  multi::array<double, 1> A(multi::index_extensions<1>{3}); assert( size(A)==3 and A[0]==0 );
}{	multi::array<double, 1> A({3}); assert( size(A)==1 and A[0]==3. );  // uses init_list
}{  multi::array<double, 1> A({3l}); assert( size(A)==1 and A[0]==3. ); // uses init_list
}{  multi::array<double, 1> A(multi::index_extensions<1>{{0, 3}}); assert( size(A)==3 and A[0]==0 );
}{  multi::array<double, 1> A({0l, 3l}); assert( size(A)==2 and A[1]==3. ); //uses init_list
}{  multi::array<double, 1, std::allocator<double>> A(multi::index_extensions<1>{2}, std::allocator<double>{}); assert( size(A)==2 );
}{  multi::array<double, 1, std::allocator<double>> A(multi::index_extensions<1>{{0, 3}}, std::allocator<double>{}); assert( size(A)==3 );
}{  multi::array<double, 2> A({2, 3}); assert( num_elements(A)==6 );
}{  multi::array<double, 2> A(multi::index_extensions<2>{{0,2}, {0,3}}); assert( num_elements(A)==6 );
}{  multi::array<double, 2, std::allocator<double>> A({2, 3}, std::allocator<double>{}); assert( num_elements(A)==6 );
}{  multi::array<double, 3> A({2, 3, 4}); assert( num_elements(A)==24 and A[1][2][3]==0 );
}{  multi::array<double, 3> A({{0, 2}, {0, 3}, {0, 4}}); assert( num_elements(A)==24 and A[1][2][3]==0 );
}{  multi::array<double, 3, std::allocator<double>> A({2, 3, 4}, std::allocator<double>{}); assert( num_elements(A)==24 );
}{  multi::array<double, 3, std::allocator<double>> A({{0, 2}, {0, 3}, {0, 4}}, std::allocator<double>{}); assert( num_elements(A)==24 );
}

 {  multi::array<double, 1> A(multi::iextensions<1>{3}, 3.1); assert( size(A)==3 and A[1]==3.1 );
}{//multi::array<double, 1> A({3}, 3.1); assert( size(A)==3 and A[1]==3.1 ); // warning in clang
}{//multi::array<double, 1> A({3l}, 3.1); assert( size(A)==3 and A[1]==3.1 ); // warning in clang
}{  multi::array<double, 1> A({0,3}, 3.1); assert( size(A)==3 and A[1]==3.1 );
}{  multi::array<double, 1> A({{0,3}}, 3.1); assert( size(A)==3 and A[1]==3.1 );
}{  multi::array<double, 1> A(multi::iextension(3), 3.1); assert( size(A)==3 and A[1]==3.1 );
}{  multi::array<double, 1> A(3l, 3.1); assert( size(A)==3 and A[1]==3.1 );
}{  multi::array<double, 1> A(3, 3.1); assert( size(A)==3 and A[1]==3.1 );
}{  multi::array<double, 1> A({0, 3}, 3.1); assert( size(A)==3 and A[1]==3.1 );
}{  multi::array<double, 2> A({2, 3}, 3.1); assert( num_elements(A)==6 and A[1][2]==3.1 );
}{  multi::array<double, 2> A({{0,2}, {0,3}}, 3.1); assert( num_elements(A)==6 and A[1][2]==3.1 );
}{  multi::array<double, 3> A({2, 3, 4}, 3.1); assert( num_elements(A)==24 and A[1][2][3]==3.1 );
}

{  
	multi::array<double, 1> A1(multi::iextension{2}, 3.1); assert( num_elements(A1)==2 and A1[1]==3.1 );
	multi::array<double, 2> A2(3, A1); assert( num_elements(A2)==6 and A2[1][2]==3.1 );
	multi::array<double, 3> A3(4, A2); assert( num_elements(A3)==24 and A3[1][2][3]==3.1 );
	
	multi::array<double, 3> B3(4, multi::array<double, 2>(3, multi::array<double, 1>({0, 2}, 3.1)));
	assert( num_elements(B3)==24 and B3[1][2][3]==3.1);
}

{  
 	std::vector<double> v1 = {1.,2.,3.};
 	multi::array<double, 1> A1(begin(v1), end(v1)); assert( size(A1)==3 and A1[2]==3. );
 	std::vector<multi::array<double, 1>> v2 = {A1, A1, A1, A1};
 	multi::array<double, 2> A2(begin(v2), end(v2)); assert( num_elements(A2)==12 and A2[1][2] == 3. );
 	std::vector<multi::array<double, 2>> v3 = {A2, A2, A2, A2, A2};
 	multi::array<double, 3> A3(begin(v3), end(v3)); assert( num_elements(A3)==60 and A3[2][1][2] == 3. );
}

 {	multi::array<double, 1> A(multi::iextensions<1>{2}); multi::array<double, 1> B=A; assert(A == B);
}{	multi::array<double, 2> A({2, 3}); multi::array<double, 2> B=A; assert(A == B);
}{	multi::array<double, 3> A({2, 3, 4}); multi::array<double, 3> B=A; assert(A == B);
}

 {	multi::array<double, 1> A(multi::iextensions<1>{2}); multi::array<double, 1> B(A, std::allocator<double>{}); assert(A == B);
}{	multi::array<double, 2> A({2, 3}); multi::array<double, 2> B(A, std::allocator<double>{}); assert(A == B);
}{	multi::array<double, 3> A({2, 3, 4}); multi::array<double, 3> B(A, std::allocator<double>{}); assert(A == B);
}

 {	multi::array<double, 1> A(multi::iextensions<1>{2}); multi::array<double, 1> B=std::move(A); assert(num_elements(B)==2 and empty(A));
}{	multi::array<double, 2> A({2, 3}); multi::array<double, 2> B=std::move(A); assert(num_elements(B)==6 and empty(A));
}{	multi::array<double, 3> A({2, 3, 4}); multi::array<double, 3> B=std::move(A); assert(num_elements(B)==24 and empty(A));
}

 {	multi::array<double, 1> A(multi::iextensions<1>{2}, 3.1); multi::array<double, 1> B(std::move(A), std::allocator<double>{}); assert(num_elements(B)==2 and B[1]==3.1 and empty(A));
}{	multi::array<double, 2> A({2, 3}, 3.1); multi::array<double, 2> B(std::move(A), std::allocator<double>{}); assert(num_elements(B)==6 and B[1][2]==3.1  and empty(A));
}{	multi::array<double, 3> A({2, 3, 4}, 3.1); multi::array<double, 3> B(std::move(A), std::allocator<double>{}); assert(num_elements(B)==24 and B[1][2][3]==3.1 and empty(A));
}

{	
	multi::array<double, 1> A1 = {1.,2.,3.}; assert(num_elements(A1)==3 and A1[1]==2.);
	multi::array<double, 2> A2 = {A1, A1, A1, A1}; assert(num_elements(A2)==12 and A2[2][1]==2.);
	multi::array<double, 3> A3 = {A2, A2, A2, A2, A2}; assert(num_elements(A3)==60 and A3[3][2][1]==2.);
}

{	
	multi::array<double, 1> A1 = {1., 2., 3.}; assert(dimensionality(A1)==1); assert(num_elements(A1)==3 and A1[1]==2.);
	multi::array<double, 1> B1 = {3., 4.}; assert(num_elements(B1)==2 and B1[1]==4.);
	multi::array<double, 1> C1 = {0, 4}; assert(num_elements(C1)==2 and C1[1]==4.);
	multi::array<double, 1> D1 = {0l, 4l}; assert(num_elements(D1)==2 and D1[1]==4.);
	multi::array<double, 1> E1 = {{0, 4}}; assert(num_elements(E1)==2 and E1[1]==4.);
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
 {	multi::array A = {1., 2., 3.}; assert( dimensionality(A)==1 and num_elements(A)==3 and A[1]==2. ); static_assert( dimensionality(A)==1 and std::rank<decltype(A)>{}==1 );
}{	multi::array A = {1., 2.};     assert( dimensionality(A)==1 and num_elements(A)==2 and A[1]==2. ); static_assert( dimensionality(A)==1 and std::rank<decltype(A)>{}==1 );
}{	multi::array A = {0, 2};       assert( dimensionality(A)==1 and num_elements(A)==2 and A[1]==2. ); static_assert( dimensionality(A)==1 and std::rank<decltype(A)>{}==1 );
}{	multi::array A = {9.};         assert( dimensionality(A)==1 and num_elements(A)==1 and A[0]==9. ); static_assert( dimensionality(A)==1 and std::rank<decltype(A)>{}==1 );
}{	multi::array A = {9};          assert( dimensionality(A)==1 and num_elements(A)==1 and A[0]==9. ); static_assert( dimensionality(A)==1 and std::rank<decltype(A)>{}==1 );
}{	multi::array A = {
		{1., 2., 3.}, 
		{4., 5., 6.}
	}; assert( dimensionality(A)==2 and num_elements(A)==6 );
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

{
    std::unique_ptr<double[]> uparr = std::make_unique<double[]>(2*3*4);
    uparr[2] = 2.;

 {	multi::array_ref<double, 1> A({24}, uparr.get()); assert( num_elements(A)==24 and A[2]==2. );
}{	multi::array_ref<double, 1> A({0}, nullptr); assert( empty(A) and num_elements(A)==0 ); 
}{	multi::array_ref<double, 1> A({0}, uparr.get()); assert( empty(A) and num_elements(A)==0 ); 
}{	multi::array_ref<double, 1> A({}, uparr.get()); assert( empty(A) and num_elements(A)==0 ); 
}{	multi::array_ref<double, 2> A({6, 4}, uparr.get()); assert( num_elements(A)==24 and A[0][2]==2. );
}{	multi::array_ref<double, 2> A({{0, 6}, {0, 4}}, uparr.get()); assert( num_elements(A)==24 and A[0][2]==2. );
}{	multi::array_ref<double, 2> A({0, 0}, uparr.get()); assert( empty(A) and num_elements(A)==0 );
}{	multi::array_ref<double, 2> A({}, uparr.get()); assert( empty(A) and num_elements(A)==0 );
}

{
	std::ptrdiff_t NX = 123, NY = 456, NZ = 789;
	std::vector<double> v(NX*NY*NZ); iota(begin(v), end(v), 0.);
	multi::array_cref<double, 3> V({NX, NY, NZ}, v.data());
	assert( V.extension(0) == NX and V.extension(1) == NY and V.extension(2) == NZ );
	assert( std::get<0>(extensions(V)) == NX and std::get<1>(extensions(V)) == NY and std::get<2>(extensions(V)) == NZ );
	for(auto i : std::get<0>(extensions(V))) assert(not V[i].empty());
}

#if __cpp_deduction_guides
 {	multi::array_ref A({24}, uparr.get()); assert( dimensionality(A)==1 and num_elements(A)==24 and A[2]==2. );
}{//multi::array_ref A({0}, nullptr); assert( dimensionality(A)==1 and empty(A) and num_elements(A)==0 ); 
}{	multi::array_ref A({0}, static_cast<double*>(nullptr)); assert( dimensionality(A)==1 and empty(A) and num_elements(A)==0 ); 
}{	multi::array_ref A({0}, uparr.get()); assert( dimensionality(A)==1 and empty(A) and num_elements(A)==0 ); 
}{//multi::array_ref A({}, uparr.get()); assert( empty(A) and num_elements(A)==0 ); 
}{	multi::array_ref A({6, 4}, uparr.get()); assert( dimensionality(A)==2 and num_elements(A)==24 and A[0][2]==2. );
}{	multi::array_ref A({{0, 6}, {0, 4}}, uparr.get()); assert( dimensionality(A)==2 and num_elements(A)==24 and A[0][2]==2. );
}{	multi::array_ref A({0, 0}, uparr.get()); assert( dimensionality(A)==2 and empty(A) and num_elements(A)==0 );
}{//multi::array_ref A({}, uparr.get()); assert( dimensionality(A)==2 and empty(A) and num_elements(A)==0 );
}
#endif

	std::vector<double> v(2*3*4);
	v[2] = 99.;
	
 {	multi::array_ref<double, 1> A({24}, v.data()); assert( num_elements(A)==24 and A[2]==99. );
}{	multi::array_ref<double, 2> A({6, 4}, v.data()); assert( num_elements(A)==24 and A[0][2]==99. );
}

 {	multi::array_ref<double, 1, std::vector<double>::iterator> A({24}, v.begin()); assert( num_elements(A)==24 and A[2]==99. );
}{	multi::array_ref<double, 2, std::vector<double>::iterator> A({6, 4}, v.begin()); assert( num_elements(A)==24 and A[0][2]==99. );
}

	std::vector<double> const vc(2*3*4, 99.);
	
 {	multi::array_ref<double, 1, double const*> A({24}, vc.data()); assert( num_elements(A)==24 and A[2]==99. );
}{	multi::array_ref<double, 2, double const*> A({6, 4}, vc.data()); assert( num_elements(A)==24 and A[0][2]==99. );
}{	multi::array_cref<double, 1> A({24}, vc.data()); assert( num_elements(A)==24 and A[2]==99. );
}{	multi::array_cref<double, 2> A({6, 4}, vc.data()); assert( num_elements(A)==24 and A[0][2]==99. );
}

	std::deque<double> q = {0.,1.,2.,3.};
 {	multi::array_ref<double, 1, std::deque<double>::iterator> A({24}, q.begin()); assert( num_elements(A)==24 and A[2]==2. );
}{	multi::array_ref<double, 2, std::deque<double>::iterator> A({6, 4}, q.begin()); assert( num_elements(A)==24 and A[0][2]==2. );
}

#if __cpp_deduction_guides
 {	multi::array_ref A({24}, v.data()); assert( dimensionality(A)==1 and num_elements(A)==24 and A[2]==99. );
}{	multi::array_ref A({6, 4}, v.data()); assert( dimensionality(A)==2 and num_elements(A)==24 and A[0][2]==99. );
}

 {	multi::array_ref A({24}, vc.data()); assert( dimensionality(A)==1 and num_elements(A)==24 and A[2]==99. );
}{	multi::array_ref A({6, 4}, vc.data()); assert( dimensionality(A)==2 and num_elements(A)==24 and A[0][2]==99. );
}

 {	multi::array_ref A({24}, v.begin()); assert( dimensionality(A)==1 and num_elements(A)==24 and A[2]==99. );
}{	multi::array_ref A({6, 4}, v.begin()); assert( dimensionality(A)==2 and num_elements(A)==24 and A[0][2]==99. );
}

 {	multi::array_ref A({24}, q.begin()); assert( dimensionality(A)==1 and num_elements(A)==24 and A[2]==2. );
}{	multi::array_ref A({6, 4}, q.begin()); assert( dimensionality(A)==2 and num_elements(A)==24 and A[0][2]==2. );
}
#endif

}

}

