#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
$CXX $0 -o $0x -lboost_unit_test_framework -lcudart -lboost_timer&&$0x&&rm $0x;exit
#endif

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi CUDA adaptor real to complex"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>
#include <boost/timer/timer.hpp>

#include "../../../../complex.hpp"
#include "../../../../array.hpp"
#include "../../../../adaptors/cuda.hpp"

namespace multi = boost::multi;
namespace cuda = multi::cuda;

BOOST_AUTO_TEST_CASE(cuda_adaptor_real_to_complex_copy){
	
	namespace mem = 
	//	multi::cuda
		multi::cuda::managed
	;

	mem::array<double              , 1> R(1<<26, 1.);
	std::cout << "elements (M) " << size(R)/1e6 << " total memory " << (size(R)*sizeof(double) + size(R)*sizeof(std::complex<double>))/1e6 << "MB" << std::endl;
	// elements (M) 67.1089 total memory 1610.61MB
	{
		boost::timer::auto_cpu_timer t; //  0.144042s wall, 0.140000s user + 0.000000s system = 0.140000s CPU (97.2%)
		mem::array<std::complex<double>, 1> C = R;
		BOOST_REQUIRE( static_cast<std::complex<double>>(C[13]) == 1. );
	}


	mem::array<std::complex<double>, 1> C = R;	

CUDA_SLOW(
	R[13] = 3.;
)

	{
		boost::timer::auto_cpu_timer t; //  0.144042s wall, 0.140000s user + 0.000000s system = 0.140000s CPU (97.2%)
		C() = R();
		BOOST_REQUIRE( static_cast<double>(R[13]) == 3. );
		BOOST_REQUIRE( static_cast<std::complex<double>>(C[13]) == 3. );
	}
	{
		boost::timer::auto_cpu_timer t; //  0.196395s wall, 0.190000s user + 0.010000s system = 0.200000s CPU (101.8%)
		C = R;
	//	BOOST_REQUIRE( static_cast<double>(R[13]) == 3. );
	//	BOOST_REQUIRE( static_cast<std::complex<double>>(C[13]) == 3. );
	}
	{
		boost::timer::auto_cpu_timer t; //  0.196395s wall, 0.190000s user + 0.010000s system = 0.200000s CPU (101.8%)
		C = R;
		BOOST_REQUIRE( static_cast<double>(R[13]) == 3. );
		BOOST_REQUIRE( static_cast<std::complex<double>>(C[13]) == 3. );
	}

	multi::array<double              , 1> Rcpu(1<<26, 1.);
	multi::array<std::complex<double>, 1> Ccpu(1<<26, 0.);
	Rcpu[13] = 3.;
	{
		boost::timer::auto_cpu_timer t; //  0.196395s wall, 0.190000s user + 0.010000s system = 0.200000s CPU (101.8%)
		Ccpu = Rcpu;
		BOOST_REQUIRE( static_cast<double>(Rcpu[13]) == 3. );
		BOOST_REQUIRE( static_cast<std::complex<double>>(Ccpu[13]) == 3. );
	}

}


