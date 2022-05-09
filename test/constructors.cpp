// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2019-2021 Alfredo A. Correa

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi constructors"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include "multi/array.hpp"

#include<complex>

namespace multi = boost::multi;

using complex = std::complex<double>;

struct multiplies_bind1st{
	explicit multiplies_bind1st(multi::array<complex, 2>&& m) : m_(std::move(m)) {} // this produces a bug in nvcc11.0
 private:
	multi::array<complex, 2> m_;
};

BOOST_AUTO_TEST_CASE(multi_construct_1d) {
	multi::static_array<double, 1> A(multi::extensions_t<1>{multi::iextension{10}}, 1.);
//	multi::static_array<double, 1> A(multi::array<double, 1>::extensions_type{10}, 1.);
	BOOST_REQUIRE( size(A) == 10 );
	BOOST_REQUIRE( A[1] == 1. );
}

BOOST_AUTO_TEST_CASE(multi_constructors_inqnvcc_bug) {
	multi::array<complex, 2> m({10, 10});
	multiplies_bind1st(std::move(m));
}

BOOST_AUTO_TEST_CASE(multi_constructors_1d) {
	{
		multi::array<double, 1> A(multi::extensions_t<1>{multi::iextension{10}});
		BOOST_REQUIRE( size(A)==10 );
	}
	{
		multi::array<double, 1> A(multi::extensions_t<1>{multi::iextension{10}}, double{});
		BOOST_REQUIRE( size(A)==10 );
		BOOST_REQUIRE( A[5]== double{} );
	}
	{
		multi::array<double, 1> A(multi::extensions_t<1>{multi::iextension{10}}, double{});
		BOOST_REQUIRE( size(A)==10 );
		BOOST_REQUIRE( A[5]== double{} );
	}
	#if defined(__cpp_deduction_guides) and not defined(__NVCC__) and not defined(__circle_build__)  // circle 170 crashes
	{
		multi::array A(multi::extensions_t<1>{{0, 10}}, double{});
		BOOST_REQUIRE( size(A)==10 );
		BOOST_REQUIRE( A[5]== double{} );
	}
	{
		multi::array A({{0, 10}}, double{});
		BOOST_REQUIRE( size(A)==10 );
		BOOST_REQUIRE( A[5]== double{} );
	}
	{
		multi::array A({10}, double{});
		BOOST_REQUIRE( size(A)==10 );
		BOOST_REQUIRE( A[5]== double{} );
	}
	{
		multi::array A(10, double{});
		BOOST_REQUIRE( size(A)==10 );
		BOOST_REQUIRE( A[5]== double{} );
	}
	#endif
}

BOOST_AUTO_TEST_CASE(multi_constructors_2d_ctad) {
#if defined(__cpp_deduction_guides) and not defined(__NVCC__) and not defined(__circle_build__)  // circle 170 crashes
	multi::array A({10, 20}, double{});
	BOOST_REQUIRE( size(A)==10 );
	BOOST_REQUIRE( A[5][6] == double{} );
#endif
}

BOOST_AUTO_TEST_CASE(multi_constructors) {
{//multi::array<double, 1> A({10}); assert(size(A)==1); // warning in clang
}{//multi::array<double, 1> A({10}, double{}); assert(size(A)==10); // warning in clang
}{//multi::array<double, 1> A({10}, double{}); assert(size(A)==10); // warning in clang
}{//multi::array<double, 1> A({10}, 0.); assert(size(A)==10); // warning in clang
}{//multi::array<double, 1> A({10}, {}); assert(size(A)==10); // error ambiguous
}{ multi::array<std::size_t, 1> A = {10}    ; BOOST_REQUIRE( size(A)==1 and A[0]==10 );
}{ multi::array<int        , 1> A = {10}    ; BOOST_REQUIRE( size(A)==1 and A[0]==10 );
}{ multi::array<double     , 1> A = {10}    ; BOOST_REQUIRE( size(A)==1 and A[0]==10 );
}{ multi::array<std::size_t, 1> A({10})     ; BOOST_REQUIRE( size(A)==1 and A[0]==10 );
}{ multi::array<int        , 1> A({10})     ; BOOST_REQUIRE( size(A)==1 and A[0]==10 );
}{ multi::array<double     , 1> A({10})     ; BOOST_REQUIRE( size(A)==1 and A[0]==10 );
//}{ multi::array<std::size_t, 1> A({{10}})   ; assert( size(A)==1 and A[0]==10 );  // clang warns about double bracked
//}{ multi::array<int        , 1> A({{10}})   ; assert( size(A)==1 and A[0]==10 );  // clang warns about double bracked
//}{ multi::array<double     , 1> A({{10}})   ; assert( size(A)==1 and A[0]==10 );  // clang warns about double bracked
}{ multi::array<std::size_t, 1> A({0, 10})  ; BOOST_REQUIRE( size(A)==2 );
}{ multi::array<int        , 1> A({0, 10})  ; BOOST_REQUIRE( size(A)==2 );
} { multi::array<double     , 1> A({0, 10})  ; BOOST_REQUIRE( size(A)==2 );
}
}

