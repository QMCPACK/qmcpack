#define BOOST_TEST_MODULE "C++ Unit Tests for Multi CUDA thrust"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include <boost/timer/timer.hpp>

#include<thrust/memory.h>

#include "../../../adaptors/thrust.hpp"
//#include "../../../memory/adaptors/cuda/cached/allocator.hpp"

#include <thrust/uninitialized_copy.h>
#include <thrust/complex.h>
#include <thrust/device_ptr.h>
#include <thrust/system/cuda/experimental/pinned_allocator.h>

namespace multi = boost::multi;

BOOST_AUTO_TEST_CASE(thrust_array) {
	multi::thrust::cuda::array<double, 2> C({2, 3});

	C[0][0] = 0. ;
	C[1][1] = 11.;
	BOOST_REQUIRE( C[1][1] == 11. );
}

BOOST_AUTO_TEST_CASE(issue_118) {

	using Allocator = thrust::device_allocator<double>;
	multi::array<double, 2, Allocator> M({3, 3}, Allocator{});

	M[1][2] = 12.;

	BOOST_REQUIRE( M[1][2] == 12. );

}

BOOST_AUTO_TEST_CASE(thrust_cpugpu_issue123){
	using T = char;

	multi::array<T, 2, thrust::cuda::allocator<T>> Devc({10240, 10240}, 'a');
	multi::array<T, 2>                             Host({10240, 10240}, 'z');

	BOOST_TEST_REQUIRE( Devc[0][0] == 'a' );
	{
		boost::timer::auto_cpu_timer t;                        // 0.010537s
		Devc = Host;
	}
	{
		boost::timer::auto_cpu_timer t;
		Devc.sliced(0, 5120) = Host.sliced(0, 5120);           //  0.005292s
	}
	{
		Host[123][456] = 'x';
		boost::timer::auto_cpu_timer t;
		Devc({0, 5120},{0, 5120}) = Host({0, 5120},{0, 5120});  // 0.002859s
	}
	BOOST_TEST_REQUIRE( Devc[123][456] == 'x' );
}

BOOST_AUTO_TEST_CASE(thrust_cpugpu_issue123_1D){
	using T = char;

	multi::array<T, 1, thrust::cuda::allocator<T>> Devc({10240*10240}, 'a');
	multi::array<T, 1>                             Host({10240*10240}, 'z');

	BOOST_TEST_REQUIRE( Devc[0] == 'a' );
	{
		boost::timer::auto_cpu_timer t;
		Devc() = Host();
	}
	{
		boost::timer::auto_cpu_timer t;
		Devc({0, 10240*10240/2}) = Host({0, 10240*10240/2});
	}
}

BOOST_AUTO_TEST_CASE(thrust_cpugpu_issue123_3D){
	using T = char;

	multi::array<T, 3, thrust::cuda::allocator<T>> Devc({1024, 1024, 1024}, 'a');
	multi::array<T, 3>                             Host({1024, 1024, 1024}, 'z');

	BOOST_TEST_REQUIRE( Devc[0][0][0] == 'a' );
	{
		boost::timer::auto_cpu_timer t;
		Devc = Host;
	}
	{
		boost::timer::auto_cpu_timer t;
		Devc({0, 512}, {0, 512}, {0, 512}) = Host({0, 512}, {0, 512}, {0, 512});
	}
}


#if 0
namespace inq{
	using complex = thrust::complex<double>;
}

BOOST_AUTO_TEST_CASE(thrust_complex_cached_1D){
	using T = inq::complex;
	multi::array<T, 1, boost::multi::memory::cuda::cached::allocator<T> > aa(10, T{1., 1.});
	multi::array<T, 1, boost::multi::memory::cuda::cached::allocator<T> > bb(10, T{2., 2.});

	bb = aa;

	BOOST_REQUIRE(( bb[0] == T{1., 1.} ));
}

BOOST_AUTO_TEST_CASE(thrust_complex_cached_without_values_1D){
	using T = inq::complex;
	multi::array<T, 1, boost::multi::memory::cuda::cached::allocator<T> > aa(10);
	multi::array<T, 1, boost::multi::memory::cuda::cached::allocator<T> > bb(10);
	BOOST_REQUIRE( aa.size() == 10 );
	BOOST_REQUIRE( bb.size() == 10 );

	bb = aa;

	BOOST_REQUIRE(( bb[0] == aa[0] ));
}

BOOST_AUTO_TEST_CASE(thrust_complex_cached_2D){
	using T = inq::complex;
	multi::array<T, 2, boost::multi::memory::cuda::cached::allocator<T> > aa({10, 20}, T{1., 1.});
	multi::array<T, 2, boost::multi::memory::cuda::cached::allocator<T> > bb({10, 20}, T{2., 2.});

	bb = aa;

	BOOST_REQUIRE(( bb[0][0] == T{1., 1.} ));
}

BOOST_AUTO_TEST_CASE(thrust_complex_cached_without_values_2D){
	using T = inq::complex;
	multi::array<T, 2, boost::multi::memory::cuda::cached::allocator<T> > aa({10, 20});
	multi::array<T, 2, boost::multi::memory::cuda::cached::allocator<T> > bb({10, 20});
	BOOST_REQUIRE( aa.size() == 10 );
	BOOST_REQUIRE( bb.size() == 10 );

	bb = aa;

	BOOST_REQUIRE(( bb[0][0] == aa[0][0] ));
}

BOOST_AUTO_TEST_CASE(array){

//{
//	multi::thrust::cuda::array<double, 2> C({2, 3});

//	C[0][0] = 0. ;
//	C[1][1] = 11.;
//	BOOST_TEST_REQUIRE( C[1][1] == 11. );
//}

//{
//	multi::array<double, 2> const H = {
//		{00., 01., 02.},
//		{10., 11., 12.},
//	};

//	BOOST_TEST_REQUIRE( H[1][1] == 11. );

//	{
//		multi::thrust::cuda::array<double, 2> C(H.extensions());
//		BOOST_REQUIRE( C.num_elements() == H.num_elements() );

//		thrust::copy_n(H.data_elements(), H.num_elements(), C.data_elements());
//		BOOST_TEST_REQUIRE( C[1][1] == 11. );
//		BOOST_REQUIRE( C == H );
//	}
//	{
//		multi::thrust::cuda::array<double, 2> C(H.extensions());
//		BOOST_REQUIRE( C.num_elements() == H.num_elements() );

//		std::copy_n(H.data_elements(), H.num_elements(), C.data_elements());
//		BOOST_TEST_REQUIRE( C[1][1] == 11. );
//		BOOST_REQUIRE( C == H );
//	}
//	{
//		multi::thrust::cuda::array<double, 2> C(H.extensions());
//		BOOST_REQUIRE( C.num_elements() == H.num_elements() );

//		std::uninitialized_copy_n(H.data_elements(), H.num_elements(), C.data_elements());
//		BOOST_TEST_REQUIRE( C[1][1] == 11. );
//		BOOST_REQUIRE( C == H );
//	}
//	{
//		multi::thrust::cuda::array<double, 2> C(H.extensions());
//		BOOST_REQUIRE( C.num_elements() == H.num_elements() );

//		thrust::uninitialized_copy_n(H.data_elements(), H.num_elements(), C.data_elements());
//		BOOST_TEST_REQUIRE( C[1][1] == 11. );
//		BOOST_REQUIRE( C == H );
//	}
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
//}

}
#endif
