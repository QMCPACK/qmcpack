#ifdef COMPILATION_INSTRUCTIONS
nvcc -x cu -O3 $0 -o $0x -lboost_unit_test_framework -D_DISABLE_CUDA_SLOW &&$0x&&rm $0x; exit
#endif

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi CUDA allocators"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include "../../cuda/allocator.hpp"
#include "../../../../array.hpp"

namespace multi = boost::multi;
namespace cuda = multi::memory::cuda;

template<class T> T what() = delete;
template<class T> T what(T&&) = delete;

BOOST_AUTO_TEST_CASE(cuda_allocators){

	multi::array<double, 1, cuda::allocator<double> > A1(200, 0.);
	using it = multi::array<double, 1, cuda::allocator<double> >::iterator;
	it sss = begin(A1);
	BOOST_REQUIRE( size(A1) == 200 );
	A1[100] = 1.;
//	what(A1.data());

	multi::array<double, 1, cuda::allocator<double>> const B1(200, 2.);
//	what(B1.data());
	BOOST_REQUIRE( B1[10] == 2. );

	A1[10] = B1[10];
	BOOST_REQUIRE( A1[10] == 2. );

	multi::array<double, 1, cuda::allocator<double>> C1(200, 0.);

//	B1[100] = 2.;
//	C1[100] = 3.;

}


