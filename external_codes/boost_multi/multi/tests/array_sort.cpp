#ifdef COMPILATION_INSTRUCTIONS
$CXX -O3 -std=c++14 -W -Wall -Wextra -Wpedantic $0 -o $0.x && $0.x $@ && rm -f $0.x; exit
#endif

#include<iostream>

#include "../array_ref.hpp"
#include "../array.hpp"

#include<algorithm> // for sort
#include<iostream> // for print
#include<cmath>
#include<vector>

namespace multi = boost::multi;
using std::cout; using std::cerr;

int main(){

	std::vector<double> v = {1.,2.,3.};
	double d2D[4][5] = {
		{150, 16, 17, 18, 19},
		{ 30,  1,  2,  3,  4}, 
		{100, 11, 12, 13, 14}, 
		{ 50,  6,  7,  8,  9} 
	};
	multi::array_ref<double, 2> d2D_ref(&d2D[0][0], {4, 5});
	
	std::stable_sort( begin(d2D_ref), end(d2D_ref) );
	std::stable_sort( d2D_ref.begin(1), d2D_ref.end(1) );

	{
		auto x = extensions(d2D_ref);
		for(auto i : std::get<0>(x)){
			for(auto j : std::get<1>(x))
				cout << d2D_ref[i][j] <<' ';
			cout <<'\n';
		}
	}

}

