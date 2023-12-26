#ifdef COMPILATION_INSTRUCTIONS
${CUDACXX:-nvcc} -std=c++17 -x cu -O3 $0 -o $0x --extended-lambda --expt-relaxed-constexpr --Werror=cross-execution-space-call -lboost_unit_test_framework -lboost_timer -Xcudafe=--display_error_number -D_DISABLE_CUDA_SLOW &&$0x&&rm $0x; exit
#endif

#define _DISABLE_CUDA_SLOW

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi CUDA array"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

//#include "../../cuda/allocator.hpp"
//#include "../../../../array.hpp"

#include "../../../../adaptors/cuda.hpp"

#include<boost/timer/timer.hpp>

namespace multi = boost::multi;
namespace cuda = multi::memory::cuda;

template<class T> T what() = delete;
template<class T> T what(T&&) = delete;

BOOST_AUTO_TEST_CASE(cuda_allocators){

	multi::array<double, 1, cuda::allocator<double> > A1(200, 0.);
	
	BOOST_REQUIRE( size(A1) == 200 );
	BEGIN_CUDA_SLOW;
	A1[100] = 1.;
	END_CUDA_SLOW;
//	what(A1.data());

	multi::array<double, 1, cuda::allocator<double>> const B1(200, 2.);
//	what(B1.data());
	BEGIN_CUDA_SLOW;
	BOOST_REQUIRE( B1[10] == 2. );
	END_CUDA_SLOW;
	
//	what(A1[10]);
//	what(B1[10]);
BEGIN_CUDA_SLOW;
	A1[10] = B1[10];
END_CUDA_SLOW;
//	BOOST_REQUIRE( A1[10] == 2. );

//	multi::array<double, 1, cuda::allocator<double>> C1(200, 0.);

//	B1[100] = 2.;
//	C1[100] = 3.;

}

BOOST_AUTO_TEST_CASE(cuda_copy_timing){
	multi::array<double, 2>::extensions_type const x = {10000, 10000};
	std::cout<<"double 2D\nsize "<< x.num_elements()*sizeof(double)/1e6 <<" MBs"<<std::endl;

	auto const cpu_times = [&]{
		multi::array<double, 2> const A(x, 999.);
		multi::array<double, 2> B(x, 888.);
		boost::timer::cpu_timer timer;
		~B() = ~A();
		BOOST_REQUIRE( B[1000][1000] == 999. );
		return timer.elapsed();
	}();
	std::cout<<"cpu "<< cpu_times.wall/1e9 <<" sec"<<std::endl;

	auto const gpu_times = [&]{
		multi::array<double, 2, cuda::allocator<double>> const A(x, 999.);
		multi::array<double, 2, cuda::allocator<double>> B(x, 888.);
		boost::timer::cpu_timer timer;
		~B() = ~A();
		BOOST_REQUIRE( B[1000][1000] == 999. );
		return timer.elapsed();
	}();
	std::cout<<"gpu "<< gpu_times.wall/1e9 <<" sec"<<std::endl;

	auto const mng_times = [&]{
		multi::array<double, 2, cuda::managed::allocator<double>> const A(x, 999.);
		multi::array<double, 2, cuda::managed::allocator<double>> B(x, 888.);
		B() = A();
		boost::timer::cpu_timer timer;
		~B() = ~A();
		BOOST_REQUIRE( B[1000][1000] == 999. );
		return timer.elapsed();
	}();
	std::cout<<"mng "<< mng_times.wall/1e9 <<" sec"<<std::endl;
}

BOOST_AUTO_TEST_CASE(cuda_copy_complex_timing){
	using complex = std::complex<double>;
	multi::array<complex, 2>::extensions_type const x = {10000, 10000};
	std::cout<<"complex 2D\nsize "<< x.num_elements()*sizeof(complex)/1e6 <<" MBs"<<std::endl;

	static_assert( std::is_trivially_copyable<std::complex<double>>{} );

	auto const cpu_times = [&]{
		multi::array<complex, 2> const A(x, 999.);
		multi::array<complex, 2> B(x, 888.);
		boost::timer::cpu_timer timer;
		B() = A();
		BOOST_REQUIRE( B[1000][1000] == 999. );
		return timer.elapsed();
	}();
	std::cout<<"cpu "<< cpu_times.wall/1e9 <<" sec"<< std::endl;

	auto const gpu_times = [&]{
		multi::array<complex, 2, cuda::allocator<complex>> const A(x, 999.);
		multi::array<complex, 2, cuda::allocator<complex>> B(x, 888.);
		boost::timer::cpu_timer timer;
		~B() = ~A();
		BOOST_REQUIRE( B[1000][1000] == 999. );
		return timer.elapsed();
	}();
	std::cout<<"gpu "<< gpu_times.wall/1e9 <<" sec"<< std::endl;

	auto const mng_times = [&]{
		multi::array<complex, 2, cuda::managed::allocator<complex>> const A(x, 999.);
		multi::array<complex, 2, cuda::managed::allocator<complex>> B(x, 888.);
		B() = A();
		boost::timer::cpu_timer timer;
		~B() = ~A();
		BOOST_REQUIRE( B[1e3][1e3] == 999. );
		return timer.elapsed();
	}();
	std::cout<<"mng "<< mng_times.wall/1e9 <<" sec"<< std::endl;
}

BOOST_AUTO_TEST_CASE(cuda_managed_empty){
	using complex = std::complex<double>;
	multi::array<complex, 2, cuda::managed::allocator<complex>> A;
	multi::array<complex, 2, cuda::managed::allocator<complex>> B = A;
	BOOST_REQUIRE( A.is_empty() );
	BOOST_REQUIRE( B.is_empty() );
	BOOST_REQUIRE( A == B );
}

BOOST_AUTO_TEST_CASE(cuda_copy_complex_timing_4d){
	using complex = std::complex<double>;
	multi::array<complex, 4>::extensions_type const x = {100, 100, 100, 100};
	std::cout<<"complex 4D\nsize "<< x.num_elements()*sizeof(complex)/1e6 <<" MBs"<<std::endl;

	static_assert( std::is_trivially_copyable<std::complex<double>>{} );

	auto const cpu_times = [&]{
		multi::array<complex, 4> const A(x, 999.);
		multi::array<complex, 4> B(x, 888.);
		boost::timer::cpu_timer timer;
		B() = A();
		BOOST_REQUIRE( B[10][10][10][10]== 999. );
		return timer.elapsed();
	}();
	std::cout<<"cpu "<< cpu_times.wall/1e9 <<" sec"<< std::endl;

	auto const gpu_times = [&]{
		multi::array<complex, 4, cuda::allocator<complex>> const A(x, 999.);
		multi::array<complex, 4, cuda::allocator<complex>> B(x, 888.);
		boost::timer::cpu_timer timer;
		~B() = ~A();
		BOOST_REQUIRE( B[10][10][10][10]== 999. );
		return timer.elapsed();
	}();
	std::cout<<"gpu "<< gpu_times.wall/1e9 <<" sec"<< std::endl;

	auto const mng_times = [&]{
		multi::array<complex, 4, cuda::managed::allocator<complex>> const A(x, 999.);
		multi::array<complex, 4, cuda::managed::allocator<complex>> B(x, 888.);
		B() = A();
		boost::timer::cpu_timer timer;
		~B() = ~A();
		BOOST_REQUIRE( B[10][10][10][10]== 999. );
		return timer.elapsed();
	}();
	std::cout<<"mng "<< mng_times.wall/1e9 <<" sec"<< std::endl;
}

