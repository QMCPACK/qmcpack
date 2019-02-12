#ifdef COMPILATION_INSTRUCTIONS
$CXX -O3 -std=c++14 -Wall -Wextra -Wpedantic -I$HOME/include $0 -o $0.x && $0.x $@ && rm -f $0.x; exit
#endif

#include<iostream>

#include "../array_ref.hpp"
#include "../array.hpp"

//#include <range/v3/algorithm/fill.hpp>
//#include <range/v3/algorithm/all_of.hpp>

#include<iostream>
#include<vector>


namespace multi = boost::multi;
using std::cout; using std::cerr;

int main(){

	multi::array<double, 2> d2D= {
		{150, 16, 17, 18, 19},
		{ 5.,  5.,  5.,  5.,  5.}, 
		{100, 11, 12, 13, 14}, 
		{ 50,  6,  7,  8,  9} 
	};
	assert( std::all_of(begin(d2D[1]), end(d2D[1]), [](auto&& e){return e == 5.;}) );

	std::fill(begin(d2D[1]), end(d2D[1]), 8.);
	assert( std::all_of(begin(d2D[1]), end(d2D[1]), [](auto&& e){return e == 8.;}) );

	std::fill(begin(d2D({0, 4}, 1)), end(d2D({0, 4}, 1)), 9.); //std::fill(begin(rotated(d2D)[1]), end(rotated(d2D)[1]), 9.);
	assert( std::all_of(begin(rotated(d2D)[1]), end(rotated(d2D)[1]), [](auto&& e){return e == 9.;}) );

//	ranges::v3::fill(d2D({0, 4}, 1), 42.);
//	ranges::v3::all_of(d2D({0, 4}, 1), [](auto&& e){return e==42.;});
}

