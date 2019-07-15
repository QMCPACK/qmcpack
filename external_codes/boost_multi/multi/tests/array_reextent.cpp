#ifdef COMPILATION_INSTRUCTIONS
$CXX catch_main.o -std=c++17 -Wall -Wextra $0 -o $0x && $0x $@ && rm $0x; exit
#endif

#include<catch2/catch.hpp>

#include "../array.hpp"

namespace multi = boost::multi;

TEST_CASE( "Array reextent", "[array]"){
	multi::array<double, 2> A({2, 3});
	REQUIRE( num_elements(A) == 6 );

	A[1][2] = 6.;
	REQUIRE( A[1][2] == 6. );

	multi::array<double, 2> C({2, 3}); 
	REQUIRE(size(C) == 2);
	REQUIRE(size(C[0]) == 3);

	A.reextent({5, 4}, 99.); 
	REQUIRE( num_elements(A)== 20 );
	REQUIRE( A[1][2] == 6. );  // reextent preserves values when it can...
	REQUIRE( A[4][3] == 99. ); // ...and gives selected value to the rest

	A = multi::array<double, 2>(extensions(A), 123.); // this is not inefficient
	REQUIRE( A[1][2] == 123. );

	clear(A); // A.clear();
	REQUIRE( num_elements(A) == 0 );
	REQUIRE( size(A) == 0 );

	A.reextent({5, 4}, 66.);
	REQUIRE( A[4][3] == 66. );		
}

