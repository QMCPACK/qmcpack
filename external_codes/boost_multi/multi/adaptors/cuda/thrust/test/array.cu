#define BOOST_TEST_MODULE "C++ Unit Tests for Multi CUDA thrust"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include "../../../../adaptors/cuda/thrust.hpp"

namespace multi = boost::multi;

template<class T> void what(T&&) = delete;

BOOST_AUTO_TEST_CASE(array){

{
	multi::thrust::cuda::array<double, 2> C({2, 3});

	C[0][0] = 0. ;
	C[1][1] = 11.;
	
	BOOST_TEST_REQUIRE( C[1][1] == 11. );
}

{
	multi::array<double, 2> const H = {
		{00., 01., 02.},
		{10., 11., 12.},
	};

	BOOST_TEST_REQUIRE( H[1][1] == 11. );

	{
		multi::thrust::cuda::array<double, 2> C(H.extensions());
		BOOST_REQUIRE( C.num_elements() == H.num_elements() );
		
		thrust::copy_n(H.data_elements(), H.num_elements(), C.data_elements());
		BOOST_TEST_REQUIRE( C[1][1] == 11. );
		BOOST_REQUIRE( C == H );
	}
	{
		multi::thrust::cuda::array<double, 2> C(H.extensions());
		BOOST_REQUIRE( C.num_elements() == H.num_elements() );
		
		std::copy_n(H.data_elements(), H.num_elements(), C.data_elements());
		BOOST_TEST_REQUIRE( C[1][1] == 11. );
		BOOST_REQUIRE( C == H );
	}
	{
		multi::thrust::cuda::array<double, 2> C(H.extensions());
		BOOST_REQUIRE( C.num_elements() == H.num_elements() );
		
		std::uninitialized_copy_n(H.data_elements(), H.num_elements(), C.data_elements());
		BOOST_TEST_REQUIRE( C[1][1] == 11. );
		BOOST_REQUIRE( C == H );
	}
	{
		multi::thrust::cuda::array<double, 2> C(H.extensions());
		BOOST_REQUIRE( C.num_elements() == H.num_elements() );
		
		what( C.data_elements() );
		thrust::uninitialized_copy_n(H.data_elements(), H.num_elements(), C.data_elements());
		BOOST_TEST_REQUIRE( C[1][1] == 11. );
		BOOST_REQUIRE( C == H );
	}
//	{
//		multi::thrust::cuda::array<double, 2> C(H.extensions());
//		BOOST_REQUIRE( C.extensions() == H.extensions() );
//		thrust::copy_n(H.begin(), H.size(), C.begin());
//		BOOST_REQUIRE( C == H );
//	}
//	{
//		multi::thrust::cuda::array<double, 2> C(H.extensions());
//		BOOST_REQUIRE( C.extensions() == H.extensions() );
//		std::copy_n(H.begin(), H.size(), C.begin());
//		BOOST_REQUIRE( C == H );
//	}
//	{
//		multi::thrust::cuda::array<double, 2> C(H.extensions());
//		C = H;
//		BOOST_REQUIRE( C == H );
//	}
//	{
//		multi::thrust::cuda::array<double, 2> C = H;
//		BOOST_REQUIRE( C == H );
//	}
}

}

