#ifdef COMPILATION_INSTRUCTIONS
$CXX -O3 -std=c++14 -Wall -Wextra -Wpedantic -Wfatal-errors $0 -o $0.x && $0.x $@ && rm -f $0.x; exit
#endif
//  (C) Copyright Alfredo A. Correa 2018.
#include "../array_ref.hpp"
#include "../array.hpp"
#include "../utility.hpp"

//#include<boost/multi_array.hpp>
#include<iostream>
#include<tuple>

using std::cout; using std::cerr;
namespace multi = boost::multi;

template<class T, typename = decltype(std::declval<T>().f())>
std::true_type has_f_aux(T const&);
std::false_type has_f_aux(...);

template<class T> struct has_f : decltype(has_f_aux(std::declval<T>())){};

struct A{
	int n;
	A() : n{5}{}
	void f() const{};
};
struct B{
	int n;
	B() : n{5}{}
};

int main(){
	assert( has_f<A>{} );
	assert( not has_f<B>{} );
	assert( not has_f<std::string>{} );
	
#if defined(__INTEL_COMPILER) or (defined(__GNUC__) and (__GNUC__ < 6))
	multi::array<double, 3> AAAA(multi::iextensions<3>{50, 50, 50});
#else
	multi::array<double, 3> AAAA({50, 50, 50});
#endif
	assert( AAAA.size() == 50 );
	assert( AAAA[0].size() == 50 );
	assert( AAAA[0][0].size() == 50 );
	assert( size(AAAA) == 50 );
	assert( size(AAAA[0]) == 50 );
	assert( size(AAAA[0][0]) == 50 );
	
	
	double DA[50][50][50];
	using multi::size;
	assert( size(DA) == 50);

	using multi::extension;
	assert(( extension(DA) == multi::index_extension{0, 50} ));
	assert(( extension(DA) == multi::iextension{0, 50} ));

	assert(( extension(DA) == multi::irange{0, 50} ));

	{
#if (defined(__GNUC__) and (__GNUC__ < 6))
	multi::array<double, 2> B(multi::iextensions<2>{50, 50});
#else
	multi::array<double, 2> B({50, 50});
#endif
	assert( size(B) == 50 );
	assert( B[0].sliced(10, 20).size() == 10 );
	assert( B(0, {10, 20}).dimensionality == 1 );
	assert( size(B(0, {10, 20})) == 10 );
	}

	multi::array<double, 2> AAA = 
		#if defined(__INTEL_COMPILER)
		(double[3][3])
		#endif
		{{1., 2., 3.}, 
		 {4., 5., 6.}, 
		 {7., 8., 9.}};
#if __GNUC__ < 6
	using extents2 = multi::iextensions<2>;
	multi::array<int, 2> A(extents2{4, 4});
#else
	multi::array<int, 2> A({4, 4});
#endif
	assert( size(A) == 4 );
	A[3][3] = 99.;
	
	decltype(A({0,2}, {0,2}))::decay_type Acopy = A({0,2}, {0,2});

	decltype(A({0,2}, {0,2})) Abb[2][2] = 
		{{A({0,2}, {0,2}), A({0, 2}, {2, 4})},
		 {A({2,4}, {0,2}), A({2, 4}, {2, 4})}};
	assert( &Abb[1][1][1][1] == &A[3][3] );

	return 0;
}

