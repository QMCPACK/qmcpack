#ifdef COMPILATION_INSTRUCTIONS
$CXX -O3 -std=c++14 -Wall -Wextra -Wpedantic `#-Wfatal-errors` $0 -o $0.x && $0.x $@ && rm -f $0.x; exit
#endif

#include<iostream>
#include<vector>
#include "../array.hpp"

namespace multi = boost::multi;

using std::cout;
void print(double const& d){cout<<d;}

template<class MultiArray> 
decltype(auto) print(MultiArray const& ma){
	cout<<"{";
	if(not ma.empty()){
		auto b = ma.begin();
		print(*b);
		std::for_each(++b, ma.end(), [](auto&& e){cout<<","; print(e);});
	}
	return cout<<"} ";
}

template<class MA>
decltype(auto) take(MA&& ma){return ma[0];}

int main(){

	{
		multi::array<double, 1> arr(100, 99.); assert(size(arr) == 100);
		assert( begin(arr) < end(arr) );

	}
	{
		multi::array<double, 2> arr({100, 100}, 99.); assert(size(arr) == 100);
		assert( cbegin(arr) < cend(arr) );
	}
	{
		std::vector<double> v(10000);
		multi::array_ref<double, 2> A(v.data(), {100, 100}); assert(size(A) == 100);
		begin(A)[4][3] = 2.; // ok 
		using multi::static_array_cast;
	//	auto const& A_const = static_array_cast<double const>(A);
	//	begin(A_const)[4][3] = 2.; // error, read only
	}
	{
		std::vector<double> dd(10000);
		multi::array_ref<double, 2, std::vector<double>::iterator> arr(begin(dd), {100, 100}); assert(size(arr) == 100);
		begin(arr)[4][3] = 2.;
	//	assert( cbegin(arr)/2 );
	//	assert( cbegin(arr) < cend(arr) );
	}
	return 0;

	multi::array<double, 3>::reverse_iterator rit;
	assert(( rit.base() == multi::array<double, 3>::reverse_iterator{}.base() ));
	assert(( multi::array<double, 3>::reverse_iterator{}.base() == multi::array<double, 3>::reverse_iterator{}.base() ));
	assert(( multi::array<double, 3>::reverse_iterator{} == multi::array<double, 3>::reverse_iterator{} ));
	assert(( multi::array<double, 3>::reverse_iterator{} == multi::array<double, 3>::reverse_iterator{} ));

	multi::array<double, 3> A =
	#if defined(__INTEL_COMPILER)
		(double[3][2][2])
	#endif
		{
			{{ 1.2,  1.1}, { 2.4, 1.}},
			{{11.2,  3.0}, {34.4, 4.}},
			{{ 1.2,  1.1}, { 2.4, 1.}}
		}
	;
	assert( begin(A) < end(A) );
	assert( cbegin(A) < cend(A) );
	assert( begin(A[0]) < end(A[0]) );

	multi::array<double, 1>::const_iterator i;
	assert( begin(A[0]) < end(A[0]) );

	assert( size(A) == 3 and size(A[0]) == 2 and size(A[0][0]) == 2 and A[0][0][1] == 1.1 );
	assert(( multi::array<double, 3>::reverse_iterator{A.begin()} == rend(A) ));

	assert( begin(A) < end(A) );
	assert( cbegin(A) < cend(A) );
//	assert( crbegin(A) < crend(A) );
//	assert( crend(A) > crbegin(A) );
	assert( end(A) - begin(A) == size(A) );
	assert( rend(A) - rbegin(A) == size(A) );

	assert( size(*begin(A)) == 2 );
	assert( size(begin(A)[1]) == 2 );
	assert( &(A[1][1].begin()[0]) == &A[1][1][0] );
	assert( &A[0][1][0] == &A[0][1][0] );
	assert( &((*A.begin())[1][0]) == &A[0][1][0] );
	assert( &((*A.begin()).operator[](1)[0]) == &A[0][1][0] );
	assert( &(A.begin()->operator[](1)[0]) == &A[0][1][0] );
	assert( &(A.begin()->operator[](1).begin()[0]) == &A[0][1][0] );
	assert( &((A.begin()+1)->operator[](1).begin()[0]) == &A[1][1][0] );
	assert( &((begin(A)+1)->operator[](1).begin()[0]) == &A[1][1][0] );
	assert( &((cbegin(A)+1)->operator[](1).begin()[0]) == &A[1][1][0] );

	multi::array<double, 3>::iterator it; assert(( it == multi::array<double, 3>::iterator{} ));
	--it;
	it = begin(A);                                    assert( it == begin(A) );
	multi::array<double, 3>::iterator it2 = begin(A); assert(it == it2);
	it = end(A);                                      assert(it != it2);
	assert(it > it2);
	multi::array<double, 3>::iterator it3{it};        assert( it3 == it );
	multi::array<double, 3>::const_iterator cit;
	cit = it3;                                        assert( cit == it3 );
	assert((begin(A) == multi::array<double, 3>::iterator{rend(A)}));
	{
		std::vector<double> vv = {1.,2.,3.};
		auto it = vv.begin();
		auto rit = vv.rend();
		assert(std::vector<double>::reverse_iterator{it} == rit);
	}
}

