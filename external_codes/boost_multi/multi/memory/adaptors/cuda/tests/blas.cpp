#ifdef COMPILATION_INSTRUCTIONS
`#nvcc -x cu --expt-relaxed-constexpr`clang++ -O3 $0 -o $0x -lboost_unit_test_framework `pkg-config --libs blas` -D_DISABLE_CUDA_SLOW -lcudart &&$0x&&rm $0x; exit
#endif

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi CUDA allocators"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include          "../../cuda/allocator.hpp"
#include "../../../../../multi/array.hpp"
#include "../../../../../multi/adaptors/blas.hpp"

namespace multi = boost::multi;
namespace cuda = multi::memory::cuda;

template<class T> void what(T&&);

BOOST_AUTO_TEST_CASE(cuda_ptr_blas){
	using complex = std::complex<double>;
	multi::array<complex, 2> A2({32, 64}); 
	A2[2][4] = complex{8., 5.};

	using multi::blas::real;
	using multi::blas::imag;

	BOOST_REQUIRE( A2[2][4] == complex(8., 5.) );
	BOOST_REQUIRE( real(A2)[2][4] == 8. );
	BOOST_REQUIRE( imag(A2)[2][4] == 5. );

// TODO: not working yet
	multi::array<complex, 2, cuda::allocator<complex>> A2_gpu = A2;
	what(multi::reinterpret_array_cast<multi::blas::Complex_<double>>(A2_gpu));
//	what(real(A2_gpu));

//	BOOST_REQUIRE( real(A2_gpu)[2][4] == 8. );
//	BOOST_REQUIRE( imag(A2_gpu)[2][4] == 5. );


	multi::array<complex, 2, cuda::managed::allocator<complex>> A2_mgpu = A2;
	BOOST_REQUIRE( std::get<0>(sizes(A2_mgpu))==32 );
	BOOST_REQUIRE( std::get<1>(sizes(A2_mgpu))==64 );

	BOOST_REQUIRE( A2_mgpu[2][4] == complex(8., 5.) );
	BOOST_REQUIRE( real(A2_mgpu)[2][4] == 8. );
	BOOST_REQUIRE( imag(A2_mgpu)[2][4] == 5. );


}

