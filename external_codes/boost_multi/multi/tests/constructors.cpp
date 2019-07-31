#ifdef COMPILATION_INSTRUCTIONS
$CXX -O3 -std=c++14 -Wall -Wextra -Wpedantic `#-Wfatal-errors` $0 -o $0.x && $0.x $@ && rm $0.x; exit
#endif

#include<iostream>
#include "../array.hpp"

#include<queue>
#include<vector>
#include<numeric> // iota

namespace multi = boost::multi;
using std::cout;

int main(){

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

// {  multi::array<double, 1> A({3}); assert( size(A)==3 and A[0]==0 );
//}
 {  multi::array<double, 1> A(multi::iextensions<1>{3}); assert( size(A)==3 and A[0]==0 );
}{  multi::array<double, 1> A(multi::array<double, 1>::extensions_type{3}); assert( size(A)==3 and A[0]==0 );
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
}{  multi::array<double, 1, std::allocator<double>> A(multi::index_extensions<1>{2}, std::allocator<double>{}); assert( size(A)==2 );
}{  multi::array<double, 1, std::allocator<double>> A(multi::index_extensions<1>{{0, 3}}, std::allocator<double>{}); assert( size(A)==3 );
}{  multi::array<double, 1, std::allocator<double>> A(multi::iextensions<1>{2}, std::allocator<double>{}); assert( size(A)==2 );
}{  multi::array<double, 1, std::allocator<double>> A(multi::iextensions<1>{{0, 3}}, std::allocator<double>{}); assert( size(A)==3 );
#if not defined(__INTEL_COMPILER) and (defined(__GNUC) and __GNU_VERSION >= 600)
}{  multi::array<double, 2> A({2, 3}); assert( num_elements(A)==6 );
#endif
}{  multi::array<double, 2> A(multi::iextensions<2>{2, 3}); assert( num_elements(A)==6 );
//}{  multi::array<double, 2> A({2, 3}); assert( num_elements(A)==6 and size(A)==2 and std::get<1>(sizes(A))==3 );
}{  multi::array<double, 2> A(multi::index_extensions<2>{{0,2}, {0,3}}); assert( num_elements(A)==6 );
#if not defined(__INTEL_COMPILER) and (defined(__GNUC__) and __GNU_VERSION__ >= 600)
}{  multi::array<double, 2, std::allocator<double>> A({2, 3}, std::allocator<double>{}); assert( num_elements(A)==6 );
#endif
}{  multi::array<double, 2, std::allocator<double>> A(multi::iextensions<2>{2, 3}, std::allocator<double>{}); assert( num_elements(A)==6 );
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
return 0;

 {  multi::array<double, 1> A(multi::iextensions<1>{3}, 3.1); assert( size(A)==3 and A[1]==3.1 );
}{//multi::array<double, 1> A({3}, 3.1); assert( size(A)==3 and A[1]==3.1 ); // warning in clang
}{//multi::array<double, 1> A({3l}, 3.1); assert( size(A)==3 and A[1]==3.1 ); // warning in clang
}{  multi::array<double, 1> A({0,3}, 3.1); assert( size(A)==3 and A[1]==3.1 );
}{  multi::array<double, 1> A(multi::iextension(3), 3.1); assert( size(A)==3 and A[1]==3.1 );
return 0;
}{  multi::array<double, 1> A(3l, 3.1); assert( size(A)==3 and A[1]==3.1 );
}{  multi::array<double, 1> A(3, 3.1); assert( size(A)==3 and A[1]==3.1 );
}{  multi::array<double, 1> A({0, 3}, 3.1); assert( size(A)==3 and A[1]==3.1 );
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

return 0;
#if 0
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
#if not __INTEL_COMPILER
}{	multi::array<double, 2> A({2, 3}); multi::array<double, 2> B=A; assert(A == B);
#endif
}{	multi::array<double, 2> A(multi::iextensions<2>{2, 3}); multi::array<double, 2> B=A; assert(A == B);
#if not __INTEL_COMPILER
}{	multi::array<double, 3> A({2, 3, 4}); multi::array<double, 3> B=A; assert(A == B);
#endif
}{	multi::array<double, 3> A(multi::iextensions<3>{2, 3, 4}); multi::array<double, 3> B=A; assert(A == B);
}

 {	multi::array<double, 1> A(multi::iextensions<1>{2}); multi::array<double, 1> B(A, std::allocator<double>{}); assert(A == B);
}{	multi::array<double, 2> A(multi::iextensions<2>{2, 3}); multi::array<double, 2> B(A, std::allocator<double>{}); assert(A == B);
#if not __INTEL_COMPILER
}{	multi::array<double, 2> A({2, 3}); multi::array<double, 2> B(A, std::allocator<double>{}); assert(A == B);
#endif
}{	multi::array<double, 3> A(multi::iextensions<3>{2, 3, 4}); multi::array<double, 3> B(A, std::allocator<double>{}); assert(A == B);
#if not __INTEL_COMPILER
}{	multi::array<double, 3> A({2, 3, 4}); multi::array<double, 3> B(A, std::allocator<double>{}); assert(A == B);
#endif
}

 {	multi::array<double, 1> A(multi::iextensions<1>{2}); multi::array<double, 1> B=std::move(A); assert(num_elements(B)==2 and empty(A));
#if not __INTEL_COMPILER
}{	multi::array<double, 2> A({2, 3}); multi::array<double, 2> B=std::move(A); assert(num_elements(B)==6 and empty(A));
#endif
}{	multi::array<double, 2> A(multi::iextensions<2>{2, 3}); multi::array<double, 2> B=std::move(A); assert(num_elements(B)==6 and empty(A));
#if not __INTEL_COMPILER
}{	multi::array<double, 3> A({2, 3, 4}); multi::array<double, 3> B=std::move(A); assert(num_elements(B)==24 and empty(A));
#endif
}{	multi::array<double, 3> A(multi::iextensions<3>{2, 3, 4}); multi::array<double, 3> B=std::move(A); assert(num_elements(B)==24 and empty(A));
}

 {	multi::array<double, 1> A(multi::iextensions<1>{2}, 3.1); multi::array<double, 1> B(std::move(A), std::allocator<double>{}); assert(num_elements(B)==2 and B[1]==3.1 and empty(A));
#if not __INTEL_COMPILER
}{	multi::array<double, 2> A({2, 3}, 3.1); multi::array<double, 2> B(std::move(A), std::allocator<double>{}); assert(num_elements(B)==6 and B[1][2]==3.1  and empty(A));
#endif
}{	multi::array<double, 2> A(multi::iextensions<2>{2, 3}, 3.1); multi::array<double, 2> B(std::move(A), std::allocator<double>{}); assert(num_elements(B)==6 and B[1][2]==3.1  and empty(A));
#if not __INTEL_COMPILER
}{	multi::array<double, 3> A({2, 3, 4}, 3.1); multi::array<double, 3> B(std::move(A), std::allocator<double>{}); assert(num_elements(B)==24 and B[1][2][3]==3.1 and empty(A));
#endif
}{	multi::array<double, 3> A(multi::iextensions<3>{2, 3, 4}, 3.1); multi::array<double, 3> B(std::move(A), std::allocator<double>{}); assert(num_elements(B)==24 and B[1][2][3]==3.1 and empty(A));
}

{	
	multi::array<double, 1> A1
		#if __INTEL_COMPILER
		= (double[3])
		#endif
		{1.,2.,3.}; 
	assert(num_elements(A1)==3 and A1[1]==2.);

	multi::array<double, 2> A2 = 
		#if __INTEL_COMPILER
		(decltype(A1)[4])
		#endif
		{A1, A1, A1, A1}; 
	assert(num_elements(A2)==12 and A2[2][1]==2.);
	multi::array<double, 3> A3 
		#if __INTEL_COMPILER
		= (decltype(A2)[5])
		#endif
		{A2, A2, A2, A2, A2}; 
	assert(num_elements(A3)==60 and A3[3][2][1]==2.);
}

{
    std::unique_ptr<double[]> uparr = std::make_unique<double[]>(2*3*4);
    uparr[2] = 2.;

 {	multi::array_ref<double, 1> A(uparr.get(), {24}); assert( num_elements(A)==24 and A[2]==2. );
}{	multi::array_ref<double, 1> A(nullptr, {0}); assert( empty(A) and num_elements(A)==0 ); 
}{	multi::array_ref<double, 1> A(uparr.get(), {0}); assert( empty(A) and num_elements(A)==0 ); 
//}{	multi::array_ref<double, 1> A(uparr.get(), {}); assert( empty(A) and num_elements(A)==0 ); // fail in icc
}{	multi::array_ref<double, 2> A(uparr.get(), {6, 4}); assert( num_elements(A)==24 and A[0][2]==2. );
}{	multi::array_ref<double, 2> A(uparr.get(), {{0, 6}, {0, 4}}); assert( num_elements(A)==24 and A[0][2]==2. );
}{	multi::array_ref<double, 2> A(uparr.get(), {0, 0}); assert( empty(A) and num_elements(A)==0 );
//}{	multi::array_ref<double, 2> A(uparr.get(), {}); assert( empty(A) and num_elements(A)==0 );
}

{
	std::ptrdiff_t NX = 123, NY = 456, NZ = 789;
	std::vector<double> v(NX*NY*NZ); iota(begin(v), end(v), 0.);
	multi::array_cref<double, 3> V(v.data(), {NX, NY, NZ});
	assert( V.extension(0) == NX and V.extension(1) == NY and V.extension(2) == NZ );
	assert( std::get<0>(extensions(V)) == NX and std::get<1>(extensions(V)) == NY and std::get<2>(extensions(V)) == NZ );
	for(auto i : std::get<0>(extensions(V))) assert(not V[i].empty());
}

#if __cpp_deduction_guides
 {	multi::array_ref A(uparr.get(), {24}); assert( multi::rank<decltype(A)>{}==1 and num_elements(A)==24 and A[2]==2. );
}{//multi::array_ref A({0}, nullptr); assert( dimensionality(A)==1 and empty(A) and num_elements(A)==0 ); 
}{	multi::array_ref A(static_cast<double*>(nullptr), {0}); assert( multi::rank<decltype(A)>{}==1 and empty(A) and num_elements(A)==0 ); 
}{	multi::array_ref A(uparr.get(), {0}); assert( multi::rank<decltype(A)>{}==1 and empty(A) and num_elements(A)==0 ); 
}{//multi::array_ref A({}, uparr.get()); assert( empty(A) and num_elements(A)==0 ); 
}{	multi::array_ref A(uparr.get(), {6, 4}); assert( multi::rank<decltype(A)>{}==2 and num_elements(A)==24 and A[0][2]==2. );
}{	multi::array_ref A(uparr.get(), {{0, 6}, {0, 4}}); assert( multi::rank<decltype(A)>{}==2 and num_elements(A)==24 and A[0][2]==2. );
}{	multi::array_ref A(uparr.get(), {0, 0}); assert( multi::rank<decltype(A)>{}==2 and empty(A) and num_elements(A)==0 );
}{//multi::array_ref A({}, uparr.get()); assert( dimensionality(A)==2 and empty(A) and num_elements(A)==0 );
}
#endif

	std::vector<double> v(2*3*4);
	v[2] = 99.;
	
 {	multi::array_ref<double, 1> A(v.data(), {24}); assert( num_elements(A)==24 and A[2]==99. );
}{	multi::array_ref<double, 2> A(v.data(), {6, 4}); assert( num_elements(A)==24 and A[0][2]==99. );
}

 {	multi::array_ref<double, 1, std::vector<double>::iterator> A(v.begin(), {24}); assert( num_elements(A)==24 and A[2]==99. );
}{	multi::array_ref<double, 2, std::vector<double>::iterator> A(v.begin(), {6, 4}); assert( num_elements(A)==24 and A[0][2]==99. );
}

	std::vector<double> const vc(2*3*4, 99.);
	
 {	multi::array_ref<double, 1, double const*> A(vc.data(), {24}); assert( num_elements(A)==24 and A[2]==99. );
}{	multi::array_ref<double, 2, double const*> A(vc.data(), {6, 4}); assert( num_elements(A)==24 and A[0][2]==99. );
}{	multi::array_cref<double, 1> A(vc.data(), {24}); assert( num_elements(A)==24 and A[2]==99. );
}{	multi::array_cref<double, 2> A(vc.data(), {6, 4}); assert( num_elements(A)==24 and A[0][2]==99. );
}

	std::deque<double> q = {0.,1.,2.,3.};
 {	multi::array_ref<double, 1, std::deque<double>::iterator> A(q.begin(), {24}); assert( num_elements(A)==24 and A[2]==2. );
}{	multi::array_ref<double, 2, std::deque<double>::iterator> A(q.begin(), {6, 4}); assert( num_elements(A)==24 and A[0][2]==2. );
}

#if __cpp_deduction_guides
 {	multi::array_ref A(v.data(), {24}); assert( multi::rank<decltype(A)>{}==1 and num_elements(A)==24 and A[2]==99. );
}{	multi::array_ref A(v.data(), {6, 4}); assert( multi::rank<decltype(A)>{}==2 and num_elements(A)==24 and A[0][2]==99. );
}

 {	multi::array_ref A(vc.data(), {24}); assert( multi::rank<decltype(A)>{}==1 and num_elements(A)==24 and A[2]==99. );
}{	multi::array_ref A(vc.data(), {6, 4}); assert( multi::rank<decltype(A)>{}==2 and num_elements(A)==24 and A[0][2]==99. );
}

 {	multi::array_ref A(v.begin(), {24}); assert( multi::rank<decltype(A)>{}==1 and num_elements(A)==24 and A[2]==99. );
}{	multi::array_ref A(v.begin(), {6, 4}); assert( multi::rank<decltype(A)>{}==2 and num_elements(A)==24 and A[0][2]==99. );
}

 {	multi::array_ref A(q.begin(), {24}); assert( multi::rank<decltype(A)>{}==1 and num_elements(A)==24 and A[2]==2. );
}{	multi::array_ref A(q.begin(), {6, 4}); assert( multi::rank<decltype(A)>{}==2 and num_elements(A)==24 and A[0][2]==2. );
}
#endif
}

cout<<"end"<<std::endl;
#endif
}

