#ifdef COMPILATION_INSTRUCTIONS//-*-indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4;-*-
#g++     -Ofast                                  -x c++  $0 -o $0x -lcudart  -lcufft `pkg-config --libs fftw3` -lboost_timer -lboost_unit_test_framework&&$0x $@&&rm $0x
#clang++ -Ofast                                  -x c++  $0 -o $0x -lcudart  -lcufft `pkg-config --libs fftw3` -lboost_timer -lboost_unit_test_framework&&$0x $@&&rm $0x
 nvcc    -O3                                     -x cu   $0 -o $0x           -lcufft `pkg-config --libs fftw3` -lboost_timer -lboost_unit_test_framework&&$0x $@&&rm $0x
 clang++ -Ofast -std=c++14 --cuda-gpu-arch=sm_60 -x cuda $0 -o $0x -lcudart  -lcufft -lfftw3_threads `pkg-config --libs fftw3` -lboost_timer -lboost_unit_test_framework&&$0x $@&&rm $0x
exit
#endif
// © Alfredo A. Correa 2020
#ifndef MULTI_ADAPTORS_FFT_HPP
#define MULTI_ADAPTORS_FFT_HPP

#include "../adaptors/fftw.hpp"
#include "../adaptors/cufft.hpp"

#if not __INCLUDE_LEVEL__
#define BOOST_TEST_MODULE "C++ Unit Tests for Multi FFT adaptor"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>
#include <boost/timer/timer.hpp>
#include <boost/config.hpp>

namespace utf = boost::unit_test;

using complex = std::complex<double>;
namespace multi = boost::multi;

using std::cout;

BOOST_AUTO_TEST_CASE(fft_combinations, *utf::tolerance(0.00001)){
	cout<< "# threads is " << multi::fftw::plan::with_nthreads() <<"\n";
	cout<<"=========================================================\n";
	cout<< BOOST_PLATFORM <<' '<< BOOST_COMPILER <<' '<< __DATE__<<'\n';

	auto const in = []{
		multi::array<complex, 4> ret({32, 90, 98, 96});
		std::generate(ret.data_elements(), ret.data_elements() + ret.num_elements(), 
			[](){return complex{std::rand()/1./RAND_MAX, std::rand()/1./RAND_MAX};}
		);
		return ret;
	}();
	std::cout<<"memory size "<< in.num_elements()*sizeof(complex)/1e6 <<" MB\n";

	multi::cuda::array<complex, 4> const in_gpu = in;
	multi::cuda::managed::array<complex, 4> const in_mng = in;

	std::vector<std::array<bool, 4>> cases = {
		{false, true , true , true }, 
		{false, true , true , false}, 
		{true , false, false, false}, 
		{true , true , false, false},
		{false, false, true , false},
		{false, false, false, false},
	};

	for(auto c : cases){
		cout<<"case: "<<std::boolalpha; 
		copy(begin(c), end(c), std::ostream_iterator<bool>{cout,", "}); cout<<"\n";

		multi::array<complex, 4> out(extensions(in));
		{
			cout<<"flops "<< multi::fftw::plan(c, in, out, multi::fft::forward).flops() <<"\n";
			boost::timer::auto_cpu_timer t{"cpu____ %ws wall, CPU (%p%)\n"};
			multi::fft::dft(c, in, out, multi::fft::forward);
		}
		{
			boost::timer::auto_cpu_timer t{"cpu_hot %ws wall, CPU (%p%)\n"};
			multi::fft::dft(c, in, out, multi::fft::forward);
		}
		multi::cuda::array<complex, 4> out_gpu(extensions(in_gpu));
		{
			boost::timer::auto_cpu_timer t{"gpu_cld %ws wall, CPU (%p%)\n"};
			multi::fft::dft(c, in_gpu   , out_gpu   , multi::fft::forward);
			BOOST_TEST( abs( static_cast<complex>(out_gpu[5][4][3][1]) - out[5][4][3][1] ) == 0. );
		}
		{
			boost::timer::auto_cpu_timer t{"gpu_hot %ws wall, CPU (%p%)\n"};
			multi::fft::dft(c, in_gpu   , out_gpu   , multi::fft::forward);
//			BOOST_TEST( abs( static_cast<complex>(out_gpu[5][4][3][1]) - out[5][4][3][1] ) == 0. );
		}
		multi::cuda::managed::array<complex, 4> out_mng(extensions(in_mng));
		{
			boost::timer::auto_cpu_timer t{"mng_cld %ws wall, CPU (%p%)\n"};
			multi::fft::dft(c, in_mng   , out_mng   , multi::fft::forward);
			BOOST_TEST( abs( out_mng[5][4][3][1] - out[5][4][3][1] ) == 0. );
		}
		{
			boost::timer::auto_cpu_timer t{"mng_hot %ws wall, CPU (%p%)\n"};
			multi::fft::dft(c, in_mng   , out_mng   , multi::fft::forward);
	//		BOOST_TEST( abs( out_mng[5][4][3][1] - out[5][4][3][1] ) == 0. );
		}
	}

}
#endif
#endif

