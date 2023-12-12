#ifdef COMPILATION_INSTRUCTIONS
(echo "#include\""$0"\"" > $0x.cpp) && clang++ `#-DNDEBUG` -O3 -std=c++14 -Wall -Wextra -Wpedantic -D_TEST_MULTI_ADAPTORS_LAPACK_TRIANGULAR -DADD_ $0x.cpp -o $0x.x -lblas && time $0x.x $@ && rm -f $0x.x $0x.cpp; exit
#endif
#ifndef MULTI_ADAPTORS_LAPACK_TRIANGULAR_HPP
#define MULTI_ADAPTORS_LAPACK_TRIANGULAR_HPP
// Alfredo A. Correa 2019 ©

#include "../../../multi/array.hpp"

namespace boost{
namespace multi{
namespace lapack{

template<class T>
struct uhermitian : public multi::array<T, 2>{
//	using multi::array<T, 2>::array;
	template<
		class MultiArray, 
		typename = decltype(multi::array<T, 2>{std::forward<MultiArray>(std::declval<MultiArray&>())})
	>
	explicit uhermitian(MultiArray&& ma) : multi::array<T, 2>{std::forward<MultiArray>(ma)}{}
	template<class Index> decltype(auto) operator[](Index i) const{
		return multi::array<T, 2>::operator[](i);
	}
	decltype(auto) operator()(index i, index j) const{
		return multi::array<T, 2>::operator[](std::min(i, j))[std::max(i, j)];
	}
};

}}}


#ifdef _TEST_MULTI_ADAPTORS_LAPACK_TRIANGULAR

#include<iostream>
#include<complex>

using std::cout;
namespace multi = boost::multi;
using complex = std::complex<double>;

int main(){
	auto const I = complex(0.,1.);

	multi::array<complex, 2> A =
		{{167.413               , 126.804 - 0.00143505*I, 125.114 - 0.1485590*I},
		 {126.804 + 0.00143505*I, 167.381               , 126.746 + 0.0327519*I},
		 {125.114 + 0.14855900*I, 126.746 - 0.0327519*I , 167.231              }}
	;
	multi::lapack::uhermitian<complex> H{std::move(A)};
	assert( &H(1, 2) == &H(2, 1) );

}
#endif
#endif

