#ifdef COMPILATION_INSTRUCTIONS
$CXX -O3 -std=c++14 -Wall -Wextra -Wpedantic $0 -o $0.x && $0.x $@ && rm -f $0.x;exit
#endif

#include "../../multi/array.hpp"

#include<iostream>
#include<vector>

namespace multi = boost::multi;
using std::cout;

multi::array_ref<double, 2> make_ref(double* p){
	return {p, {5, 7}};
}

int main(){

	{
		std::vector<double> v(5*7, 99.);

		multi::array<double, 2> A{{5, 7}, 33.};
		make_ref(v.data()) = A;

		assert( v[9] == 33. );
	}
	{
		std::vector<double> v(5*7, 99.), w(5*7, 33.);

		multi::array_ref<double, 2> B{w.data(), {5, 7}};
		make_ref(v.data()) = B;
		make_ref(v.data()) = B.sliced(0,5);

		assert( v[9] == 33. );
	}
	{
		std::vector<double> v(5*7, 99.), w(5*7, 33.);

		make_ref(v.data()) = make_ref(w.data());

		assert( v[9] == 33. );
	}
}

