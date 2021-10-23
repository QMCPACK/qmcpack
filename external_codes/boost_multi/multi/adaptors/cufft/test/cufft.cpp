#ifdef COMPILATION// -*-indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4;-*-
$CXX $0 -o $0x -lcudart -lcufft `pkg-config --libs fftw3` -lboost_unit_test_framework&&$0x&&rm $0x;exit
#endif
// © Alfredo A. Correa 2020-2021

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi cuFFT adaptor"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include <boost/timer/timer.hpp>

#include "../../../adaptors/cuda.hpp"
#include "../../../adaptors/fftw.hpp"
#include "../../../adaptors/cufft.hpp"

#include<complex>
#include<thrust/complex.h>
#include "../../../complex.hpp"

#include<cuda_runtime.h> // cudaDeviceSynchronize

#include<iostream>

namespace multi = boost::multi;
using complex = std::complex<double>;
namespace utf = boost::unit_test;


template <class T>
__attribute__((always_inline)) inline void DoNotOptimize(const T &value) {
	asm volatile("" : "+m"(const_cast<T &>(value)));
}

struct watch : private std::chrono::high_resolution_clock{
	std::string label_; time_point  start_;
	watch(std::string label ="") : label_{label}, start_{}{
		cudaDeviceSynchronize();
		start_ = now();
	}
	~watch(){
		cudaDeviceSynchronize();
		auto const count = std::chrono::duration<double>(now() - start_).count();
		std::cerr<< label_<<": "<< count <<" sec"<<std::endl;
	}
};

constexpr complex I{0, 1};

//#if 1

BOOST_AUTO_TEST_CASE(cufft_2D, *boost::unit_test::tolerance(0.0001)){

//	boost::multi::cuda::allocator<std::complex<double>>::const_void_pointer cvp1 = nullptr;
//	std::allocator_traits<boost::multi::cuda::allocator<std::complex<double>>>::const_void_pointer cvp2 = nullptr;

	multi::array<complex, 2> const in_cpu = {
		{ 1. + 2.*I, 9. - 1.*I, 2. + 4.*I},
		{ 3. + 3.*I, 7. - 4.*I, 1. + 9.*I},
		{ 4. + 1.*I, 5. + 3.*I, 2. + 4.*I},
		{ 3. - 1.*I, 8. + 7.*I, 2. + 1.*I},
		{ 31. - 1.*I, 18. + 7.*I, 2. + 10.*I}
	};
	multi::array<complex, 2> fw_cpu(extensions(in_cpu));
	multi::fftw::dft(in_cpu, fw_cpu, multi::fftw::forward);

	multi::cuda::array<complex, 2> const in_gpu = in_cpu;
	multi::cuda::array<complex, 2> fw_gpu(extensions(in_gpu));
	multi::cufft::dft(in_gpu, fw_gpu, multi::cufft::forward);

	BOOST_TEST( std::imag(static_cast<complex>(fw_gpu[3][2]) - fw_cpu[3][2]) == 0. );

	auto fw2_gpu = multi::cufft::dft(in_gpu, multi::cufft::forward);
	BOOST_TEST( std::imag(static_cast<complex>(fw2_gpu[3][1]) - fw_cpu[3][1]) == 0. );

	multi::cuda::managed::array<complex, 2> const in_mng = in_cpu;
	multi::cuda::managed::array<complex, 2> fw_mng(extensions(in_gpu));
	multi::cufft::dft(in_mng, fw_mng, multi::cufft::forward);

	BOOST_TEST( std::imag(fw_mng[3][2] - fw_cpu[3][2]) == 0. );

//	auto fw2_mng = multi::fftw::dft(in_mng, multi::fftw::forward);
//	BOOST_TEST( std::imag(fw2_mng[3][1] - fw_cpu[3][1]) == 0. );

}

BOOST_AUTO_TEST_CASE(cufft_3D_timing, *boost::unit_test::tolerance(0.0001)){

	auto x = std::make_tuple(300, 300, 300);
	{
		multi::array<complex, 3> const in_cpu(x, 10.); 
		BOOST_ASSERT( in_cpu.num_elements()*sizeof(complex) < 2e9 );
		multi::array<complex, 3> fw_cpu(extensions(in_cpu), 99.);
		{
		//	boost::timer::auto_cpu_timer t;  // 1.041691s wall, 1.030000s user + 0.000000s system = 1.030000s CPU (98.9%)
			multi::fftw::dft(in_cpu, fw_cpu, multi::fftw::forward);
			BOOST_TEST( fw_cpu[8][9][10] != 99. );
		}
	}
	{
		multi::cuda::array<complex, 3> const in_gpu(x, 10.); 
		multi::cuda::array<complex, 3> fw_gpu(extensions(in_gpu), 99.);
		{
		//	boost::timer::auto_cpu_timer t; //  0.208237s wall, 0.200000s user + 0.010000s system = 0.210000s CPU (100.8%)
			multi::cufft::dft(in_gpu, fw_gpu, multi::fftw::forward);

			BOOST_TEST( static_cast<complex>(fw_gpu[8][9][10]) != 99. );
		}
	}
	{
		multi::cuda::managed::array<complex, 3> const in_gpu(x, 10.); 
		multi::cuda::managed::array<complex, 3> fw_gpu(extensions(in_gpu), 99.);
		{
		//	boost::timer::auto_cpu_timer t; //  0.208237s wall, 0.200000s user + 0.010000s system = 0.210000s CPU (100.8%)
			multi::cufft::dft(in_gpu, fw_gpu, multi::cufft::forward);
		//	BOOST_TEST( fw_gpu[8][9][10].operator complex() != 99. );
		}
		{
		//	boost::timer::auto_cpu_timer t; //  0.208237s wall, 0.200000s user + 0.010000s system = 0.210000s CPU (100.8%)
			multi::cufft::dft(in_gpu, fw_gpu, multi::cufft::forward);
		//	BOOST_TEST( fw_gpu[8][9][10].operator complex() != 99. );
		}
	}
}

BOOST_AUTO_TEST_CASE(cufft_combinations, *utf::tolerance(0.00001)){

	auto const in = []{
		multi::array<complex, 4> ret({32, 90, 98, 96});
		std::generate(ret.data_elements(), ret.data_elements() + ret.num_elements(), 
			[](){return complex{std::rand()*1./RAND_MAX, std::rand()*1./RAND_MAX};}
		);
		return ret;
	}();
	std::clog<<"memory size "<< in.num_elements()*sizeof(complex)/1e6 <<" MB\n";

	multi::cuda::array<complex, 4> const in_gpu = in;
	multi::cuda::managed::array<complex, 4> const in_mng = in;

	using std::clog;
	for(auto c : std::vector<std::array<bool, 4>>{
		{false, true , true , true }, 
		{false, true , true , false}, 
		{true , false, false, false}, 
		{true , true , false, false},
		{false, false, true , false},
		{false, false, false, false},
	}){
		std::clog<<"case "; copy(begin(c), end(c), std::ostream_iterator<bool>{std::clog,", "}); std::clog<<std::endl;
		multi::array<complex, 4> out = in;
		multi::array<complex, 4> in_rw = in;
		[&, _ = watch{"cpu_opl "}]{
			multi::fftw::dft_forward(c, in, out);
		}();
		[&, _ = watch{"cpu_ipl "}]{
			multi::fftw::dft(c, in_rw, multi::fftw::forward);
			BOOST_TEST( abs( static_cast<multi::complex<double>>(in_rw[5][4][3][1]) - multi::complex<double>(out[5][4][3][1]) ) == 0. );			
		}();
		{
			multi::array<complex, 4> in_rw2 = in;
			[&, _ = watch{"cpu_mov "}]{
				multi::array<complex, 4> const out_mov = multi::fftw::dft_forward(c, std::move(in_rw2));
			//	what(out_mov);
				BOOST_TEST( abs( static_cast<multi::complex<double>>(out_mov[5][4][3][1]) - multi::complex<double>(out[5][4][3][1]) ) == 0. );			
				BOOST_REQUIRE( is_empty(in_rw2) );
				BOOST_REQUIRE( extensions(out_mov) == extensions(in) );
			}();
		}


		[&, _ = watch{"cpu_new "}]{
			auto const out_cpy = multi::fftw::dft_forward(c, in);
			BOOST_TEST( abs( static_cast<multi::complex<double>>(out_cpy[5][4][3][1]) - multi::complex<double>(out[5][4][3][1]) ) == 0. );
		}();
		multi::cuda::array<complex, 4> out_gpu(extensions(in_gpu));
		[&, _ = watch{"gpu_opl "}]{
			multi::cufft::dft(c, in_gpu   , out_gpu, multi::cufft::forward);
			BOOST_TEST( abs( static_cast<complex>(out_gpu[5][4][3][1]) - out[5][4][3][1] ) == 0. );
		}();
		{
			multi::cuda::array<complex, 4> in_rw_gpu = in_gpu;
			[&, _ = watch{"gpu_ipl "}]{
				multi::cufft::dft(c, in_rw_gpu, multi::cufft::forward);
				BOOST_TEST( abs( static_cast<complex>(in_rw_gpu[5][4][3][1]) - out[5][4][3][1] ) == 0. );
			}();
		}
		{
			multi::cuda::array<complex, 4> in_rw_gpu = in_gpu;
			[&, _ = watch{"gpu_mov "}]{
				multi::cuda::array<complex, 4> const out_mov = multi::cufft::dft_forward(c, std::move(in_rw_gpu));
			//	BOOST_REQUIRE( in_rw_gpu.empty() );
			//	BOOST_TEST( abs( static_cast<complex>(out_mov[5][4][3][1]) - out[5][4][3][1] ) == 0. );
			}();
		}
		{
			multi::cuda::array<complex, 4> in_rw_gpu = in_gpu;
			[&, _ = watch{"gpu_mov "}]{
				multi::cuda::array<complex, 4> out_mov = std::move(in_rw_gpu);
				multi::cufft::dft(c, out_mov, multi::cufft::forward);
			//	BOOST_REQUIRE( in_rw_gpu.empty() );
			//	BOOST_TEST( abs( static_cast<complex>(out_mov[5][4][3][1]) - out[5][4][3][1] ) == 0. );
			}();
		}
		cudaDeviceSynchronize();
		[&, _ = watch{"gpu_new "}]{
			multi::cuda::array<complex, 4> const out_cpy = multi::cufft::dft(c, in_gpu, multi::cufft::forward);
		}();
		multi::cuda::managed::array<complex, 4> out_mng(extensions(in_mng));
		[&, _ = watch{"mng_cld "}]{
			multi::cufft::dft(c, in_mng, out_mng, multi::cufft::forward);
			BOOST_TEST( abs( out_mng[5][4][3][1] - out[5][4][3][1] ) == 0. );
		}();
		[&, _ = watch{"mng_hot "}]{
			multi::cufft::dft(c, in_mng   , out_mng, multi::cufft::forward);
			BOOST_TEST( abs( out_mng[5][4][3][1] - out[5][4][3][1] ) == 0. );
		}();
		[&, _ = watch{"mng_new "}]{
			auto const out_mng = multi::cufft::dft(c, in_mng, multi::cufft::forward);
			BOOST_TEST( abs( out_mng[5][4][3][1] - out[5][4][3][1] ) == 0. );
		}();
	}
	std::clog<<"cache size "
		<< multi::cufft::plan::cache<1>().size() <<' '
		<< multi::cufft::plan::cache<2>().size() <<' '
		<< multi::cufft::plan::cache<3>().size() <<' '
		<< multi::cufft::plan::cache<4>().size() <<' '
	<<std::endl;

}

BOOST_AUTO_TEST_CASE(cufft_many_3D, *utf::tolerance(0.00001) ){

	auto const in_cpu = []{
		multi::array<complex, 4> ret({45, 18, 32, 16});
		std::generate(
			ret.data_elements(), ret.data_elements() + ret.num_elements(), 
			[](){return complex{std::rand()*1./RAND_MAX, std::rand()*1./RAND_MAX};}
		);
		return ret;
	}();

	multi::cuda::array<complex, 4> const in = in_cpu;
	multi::cuda::array<complex, 4>       out(extensions(in));

#if 0
	multi::cufft::many_dft(begin(unrotated(in)), end(unrotated(in)), begin(unrotated(out)), +1);

	multi::array<complex, 4> out_cpu(extensions(in));
	multi::fft::many_dft(begin(unrotated(in_cpu)), end(unrotated(in_cpu)), begin(unrotated(out_cpu)), +1);

	BOOST_TEST( imag( static_cast<complex>(out[5][4][3][2]) - out_cpu[5][4][3][2]) == 0. );
#endif
}

BOOST_AUTO_TEST_CASE(cufft_4D, *utf::tolerance(0.00001) ){
	auto const in = []{
		multi::array<complex, 3> ret({10, 10, 10});
		std::generate(ret.data_elements(), ret.data_elements() + ret.num_elements(), 
			[](){return complex{std::rand()*1./RAND_MAX, std::rand()*1./RAND_MAX};}
		);
		return ret;
	}();

	multi::array<complex, 3> out(extensions(in));
//	multi::fftw::dft({true, false, true}, in, out, multi::fftw::forward);
	multi::fftw::many_dft(begin(in<<1), end(in<<1), begin(out<<1), multi::fftw::forward);

	multi::cuda::array<complex, 3> in_gpu = in;
	multi::cuda::array<complex, 3> out_gpu(extensions(in));

//	multi::cufft::dft({true, false, true}, in_gpu, out_gpu, multi::fft::forward);//multi::cufft::forward);	
	multi::cufft::many_dft(begin(in_gpu<<1), end(in_gpu<<1), begin(out_gpu<<1), multi::fftw::forward);
	BOOST_TEST( imag( static_cast<complex>(out_gpu[5][4][3]) - out[5][4][3]) == 0. );	
}

