#ifdef COMPILATION_INSTRUCTIONS//-*-indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4;-*-
$CXX $0 -o $0x -lcudart -lboost_timer -lboost_unit_test_framework&&$0x&&rm $0x;exit
#endif
// © Alfredo A. Correa 2020
#define BOOST_TEST_MODULE "C++ Unit Tests for Multi cuda adaptor"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>
namespace utf = boost::unit_test;
#include <boost/timer/timer.hpp>

#include "../../../adaptors/cuda.hpp"

#include<complex>

namespace boost::multi::memory::cuda{
template<class T1, class P1, class T2, class P2>
void copy(array_iterator<T1, 1, boost::multi::memory::cuda::ptr<P1>>, array_iterator<T1, 1, boost::multi::memory::cuda::ptr<P1>>, array_iterator<T2, 1, boost::multi::memory::cuda::ptr<P2>>){
	assert(0);
}
//std::copy<boost::multi::array_iterator<double, 1, boost::multi::memory::cuda::ptr<const double, const double *>, boost::multi::memory::cuda::ref<const double> >, boost::multi::array_iterator<double, 1, boost::multi::memory::cuda::ptr<double, double *>, boost::multi::memory::cuda::ref<double> > >
}

namespace multi = boost::multi;

BOOST_AUTO_TEST_CASE(cudart_double, *utf::tolerance(0.00001)*utf::timeout(10)){

	auto const in = []{
		multi::array<double, 4> r({32, 90, 98, 96});
		std::generate(data_elements(r), data_elements(r)+num_elements(r), &std::rand);
		return r;
	}();
	std::cout<<"memory size "<< in.num_elements()*sizeof(decltype(in)::element)/1e6 <<" MB\n";

	{
		boost::timer::auto_cpu_timer t{"%ws wall, CPU (%p%)\n"};
		multi::cuda::array<double, 4> const in_gpu = in;

		multi::array<double, 4> const in_cpy = in_gpu;
		BOOST_REQUIRE( in == in_cpy );
	}
	{
		boost::timer::auto_cpu_timer t{"%ws wall, CPU (%p%)\n"};
		multi::cuda::array<double, 4> const in_gpu = in;
	}
	{
		multi::cuda::array<double, 4> const in_gpu = in;
		multi::cuda::array<double, 4> out_gpu = in;
		boost::timer::auto_cpu_timer t{"copy assign gpu____ %ws wall, CPU (%p%)\n"};
		out_gpu = in_gpu;
		auto c = static_cast<double>(out_gpu[1][2][3][4]); (void)c;

		(out_gpu << 1) = (in_gpu << 1);
	}
	{
		multi::cuda::managed::array<double, 4> const in_mng = in;
		multi::cuda::managed::array<double, 4> out_mng = in;
		{
			boost::timer::auto_cpu_timer t{"copy assign mng____ %ws wall, CPU (%p%)\n"};
			out_mng = in_mng;
			auto c = static_cast<double>(out_mng[1][2][3][4]); (void)c;
		}
		{
			boost::timer::auto_cpu_timer t{"copy assign mng_hot %ws wall, CPU (%p%)\n"};
			out_mng = in_mng;
			auto c = static_cast<double>(out_mng[1][2][3][4]); (void)c;
		}
		{
			boost::timer::auto_cpu_timer t{"copy assign mng loop %ws wall, CPU (%p%)\n"};
			out_mng() = in_mng();
			auto c = static_cast<double>(out_mng[1][2][3][4]); (void)c;
		}
	}
}

BOOST_AUTO_TEST_CASE(cudart_complex, *utf::tolerance(0.00001)*utf::timeout(10)){

	using complex = std::complex<double>;

	auto const in = []{
		multi::array<complex, 4> r({32, 90, 98, 96});
		std::generate(data_elements(r), data_elements(r)+num_elements(r), &std::rand);
		return r;
	}();
	std::cout<<"memory size "<< in.num_elements()*sizeof(decltype(in)::element)/1e6 <<" MB\n";

	{
		boost::timer::auto_cpu_timer t{"%ws wall, CPU (%p%)\n"};
		multi::cuda::array<complex, 4> const in_gpu = in;
	}
	{
		boost::timer::auto_cpu_timer t{"%ws wall, CPU (%p%)\n"};
		multi::cuda::array<complex, 4> const in_gpu = in;
	}

}
