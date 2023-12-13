// © Alfredo A. Correa 2020-2023

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi cuFFT adaptor"
#include<boost/test/unit_test.hpp>

#include <boost/timer/timer.hpp>

// #include "../../../adaptors/cuda.hpp"
#include "../../../adaptors/fft.hpp"
#include "../../../adaptors/fftw.hpp"
#include "../../../adaptors/cufft.hpp"
#include "../../../adaptors/thrust.hpp"

#include<thrust/complex.h>
#include "../../../complex.hpp"

#include<cuda_runtime.h> // cudaDeviceSynchronize

#include <chrono>
#include <complex>
#include <iostream>
#include <random>

namespace multi = boost::multi;
using complex = thrust::complex<double>;
namespace utf = boost::unit_test;
complex const I{0.0, 1.0};

template<>
constexpr bool multi::force_element_trivial_default_construction<thrust::complex<double>> = true;

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

BOOST_AUTO_TEST_CASE(cufft_2D, *boost::unit_test::tolerance(0.0001)){

	using complex = thrust::complex<double>;  // this can't be std::complex<double> in the gpu

	auto const in_cpu = multi::array<complex, 2>{
		{ 1.0 + 2.0*I, 9.0 - 1.0*I, 2.0 + 4.0*I},
		{ 3.0 + 3.0*I, 7.0 - 4.0*I, 1.0 + 9.0*I},
		{ 4.0 + 1.0*I, 5.0 + 3.0*I, 2.0 + 4.0*I},
		{ 3.0 - 1.0*I, 8.0 + 7.0*I, 2.0 + 1.0*I},
		{ 31.0 - 1.0*I, 18.0 + 7.0*I, 2.0 + 10.0*I}
	};

	{
		auto       fw_cpu = multi::array<complex, 2>(extensions(in_cpu));
		multi::fftw::dft_forward({true, true}, in_cpu, fw_cpu);

		auto const in_gpu = multi::thrust::cuda::array<complex, 2>{in_cpu};
		auto       fw_gpu = multi::thrust::cuda::array<complex, 2>(extensions(in_gpu));

		BOOST_TEST( fw_cpu[3][2].real() != 0.0 );
		BOOST_TEST( fw_cpu[3][2].imag() != 0.0 );

		multi::cufft::plan<2>({true, true}, in_gpu.layout(), fw_gpu.layout())
			.execute(in_gpu.base(), fw_gpu.base(), multi::cufft::forward);

		BOOST_TEST( (complex(fw_gpu[3][2]) - fw_cpu[3][2]).real() == 0.0 );
		BOOST_TEST( (complex(fw_gpu[3][2]) - fw_cpu[3][2]).imag() == 0.0 );
	}
	{
		auto       fw_cpu = multi::array<complex, 2>(extensions(in_cpu));
		multi::fftw::dft({false, true}, in_cpu, fw_cpu, multi::fftw::forward);

		auto const in_gpu = multi::thrust::cuda::array<complex, 2>{in_cpu};
		auto       fw_gpu = multi::thrust::cuda::array<complex, 2>(extensions(in_gpu));

		BOOST_TEST( fw_cpu[3][2].real() != 0.0 );
		BOOST_TEST( fw_cpu[3][2].imag() != 0.0 );

		multi::cufft::plan<2>({false, true}, in_gpu.layout(), fw_gpu.layout())
			.execute(in_gpu.base(), fw_gpu.base(), multi::cufft::forward);

		BOOST_TEST( (complex(fw_gpu[3][2]) - fw_cpu[3][2]).real() == 0.0 );
		BOOST_TEST( (complex(fw_gpu[3][2]) - fw_cpu[3][2]).imag() == 0.0 );
	}
	{
		auto       fw_cpu = multi::array<complex, 2>(extensions(in_cpu));
		multi::fftw::dft({false, true}, in_cpu, fw_cpu, multi::fftw::forward);

		auto const in_gpu = multi::thrust::cuda::array<complex, 2>{in_cpu};
		auto       fw_gpu = multi::thrust::cuda::array<complex, 2>(extensions(in_gpu));

		BOOST_TEST( fw_cpu[3][2].real() != 0.0 );
		BOOST_TEST( fw_cpu[3][2].imag() != 0.0 );

		for(int i = 0; i != in_gpu.size(); ++i) {
			multi::cufft::plan<1>({true}, in_gpu[i].layout(), fw_gpu[i].layout())
				.execute(in_gpu[i].base(), fw_gpu[i].base(), multi::cufft::forward);
		}

		BOOST_TEST( (complex(fw_gpu[3][2]) - fw_cpu[3][2]).real() == 0.0 );
		BOOST_TEST( (complex(fw_gpu[3][2]) - fw_cpu[3][2]).imag() == 0.0 );
	}
	{
		auto       fw_cpu = multi::array<complex, 2>(extensions(in_cpu));
		multi::fftw::dft({false, true}, in_cpu, fw_cpu, multi::fftw::forward);

		auto const in_gpu = multi::thrust::cuda::array<complex, 2>{in_cpu};
		auto       fw_gpu = multi::thrust::cuda::array<complex, 2>(extensions(in_gpu));
		auto       fw_gpu2 = multi::thrust::cuda::array<complex, 2>(extensions(in_gpu));
		auto       fw_gpu3 = multi::thrust::cuda::array<complex, 2>(extensions(in_gpu));

		BOOST_TEST( fw_cpu[3][2].real() != 0.0 );
		BOOST_TEST( fw_cpu[3][2].imag() != 0.0 );

		for(int i = 0; i != in_gpu.size(); ++i) {
			multi::cufft::plan<1>({true}, in_gpu[i].layout(), fw_gpu[i].layout())
				.execute(in_gpu[i].base(), fw_gpu[i].base(), multi::cufft::forward);
		}

		multi::cufft::plan<2>({false, true}, in_gpu.layout(), fw_gpu2.layout())
			.execute(in_gpu.base(), fw_gpu2.base(), multi::cufft::forward);

		BOOST_TEST( (complex(fw_gpu[3][2]) - fw_cpu[3][2]).real() == 0.0 );
		BOOST_TEST( (complex(fw_gpu[3][2]) - fw_cpu[3][2]).imag() == 0.0 );

		BOOST_TEST( (complex(fw_gpu[3][2]) - complex(fw_gpu2[3][2])).real() == 0.0 );
		BOOST_TEST( (complex(fw_gpu[3][2]) - complex(fw_gpu2[3][2])).imag() == 0.0 );
	}
	{
		auto       fw_cpu = multi::array<complex, 2>(extensions(in_cpu));
		multi::fftw::dft({false, true}, in_cpu, fw_cpu, multi::fftw::forward);

		auto const in_gpu = multi::thrust::cuda::array<complex, 2>{in_cpu};
		auto const fw_gpu = multi::cufft::dft({false, true}, in_gpu, multi::cufft::forward);

		BOOST_TEST( fw_cpu[3][2].real() != 0.0 );
		BOOST_TEST( fw_cpu[3][2].imag() != 0.0 );

		BOOST_TEST( (complex(fw_gpu[3][2]) - fw_cpu[3][2]).real() == 0.0 );
		BOOST_TEST( (complex(fw_gpu[3][2]) - fw_cpu[3][2]).imag() == 0.0 );

		BOOST_TEST( (complex(fw_gpu[2][3]) - fw_cpu[2][3]).real() == 0.0 );
		BOOST_TEST( (complex(fw_gpu[2][3]) - fw_cpu[2][3]).imag() == 0.0 );
	}
	{
		auto       fw_cpu = multi::array<complex, 2>(extensions(in_cpu));
		multi::fftw::dft({true, false}, in_cpu, fw_cpu, multi::fftw::forward);

		auto const in_gpu = multi::thrust::cuda::array<complex, 2>{in_cpu};
		auto const fw_gpu = multi::cufft::dft({true, false}, in_gpu, multi::cufft::forward);

		BOOST_TEST( fw_cpu[3][2].real() != 0.0 );
		BOOST_TEST( fw_cpu[3][2].imag() != 0.0 );

		BOOST_TEST( (complex(fw_gpu[3][2]) - fw_cpu[3][2]).real() == 0.0 );
		BOOST_TEST( (complex(fw_gpu[3][2]) - fw_cpu[3][2]).imag() == 0.0 );

		BOOST_TEST( (complex(fw_gpu[2][3]) - fw_cpu[2][3]).real() == 0.0 );
		BOOST_TEST( (complex(fw_gpu[2][3]) - fw_cpu[2][3]).imag() == 0.0 );
	}
}

BOOST_AUTO_TEST_CASE(check_thrust_complex_vs_std_complex, *boost::unit_test::tolerance(0.0001)){

	multi::array<std   ::complex<double>, 1> const s_in = {1.0 + I*2.0, 2.0 + I*3.0, 3.0 + I*4.0};
	multi::array<thrust::complex<double>, 1> const t_in = {1.0 + I*2.0, 2.0 + I*3.0, 3.0 + I*4.0};

	multi::array<std   ::complex<double>, 1>       s_out(s_in.extensions());
	multi::array<thrust::complex<double>, 1>       t_out(t_in.extensions());

	multi::fftw::plan::forward({true}, s_in.base(), s_in.layout(), s_out.base(), s_out.layout()).execute(s_in.base(), s_out.base());
	multi::fftw::plan::forward({true}, t_in.base(), t_in.layout(), t_out.base(), t_out.layout()).execute(t_in.base(), t_out.base());

	BOOST_REQUIRE( std::equal(s_out.begin(), s_out.end(), t_out.begin()) );
}

BOOST_AUTO_TEST_CASE(small_1D_cpu_vs_cpu, *boost::unit_test::tolerance(0.0001)){

	multi::array<thrust::complex<double>, 1> const cpu_in = {1.0 + I*2.0, 2.0 + I*3.0, 3.0 + I*4.0};
	multi::thrust::cuda::array<thrust::complex<double>, 1> const gpu_in = {1.0 + I*2.0, 2.0 + I*3.0, 3.0 + I*4.0};

	multi::array<thrust::complex<double>, 1> cpu_out(cpu_in.extensions());
	multi::thrust::cuda::array<thrust::complex<double>, 1> gpu_out(gpu_in.extensions());

	multi::fftw::plan::forward({true}, cpu_in.base(), cpu_in.layout(), cpu_out.base(), cpu_out.layout()).execute        (cpu_in.base(), cpu_out.base());
	multi::cufft::plan<1>     ({true},                gpu_in.layout(),                 gpu_out.layout()).execute_forward(gpu_in.base(), gpu_out.base());
}

BOOST_AUTO_TEST_CASE(cufft_1D_combinations, *boost::unit_test::tolerance(0.0001)){

	using complex = thrust::complex<double>;  // this can't be std::complex<double> in the gpu

	auto const in_cpu = std::invoke([]{
		multi::array<complex, 1> ret({128}, complex{});
		std::default_random_engine generator;
		std::uniform_real_distribution<double> distribution(1.0, 88.0);

		std::generate(
			reinterpret_cast<double*>(ret.data_elements()),
			reinterpret_cast<double*>(ret.data_elements() + ret.num_elements()), [&]{return distribution(generator);}
		);
		return ret;
	});

	for(auto c : std::vector<std::array<bool, 1>>{
		{true} //,
		// {false},
	}){
		auto const in_gpu = multi::thrust::cuda::array<complex, 1>{in_cpu};

		for(auto const idx : extension(in_cpu)) {
			std::cout << "A: " << idx << ": " << in_cpu[idx] << ", " << in_gpu[idx] << std::endl;
		}

		BOOST_TEST( complex(in_gpu[31]).real() == in_cpu[31].real() );
		BOOST_TEST( complex(in_gpu[31]).imag() == in_cpu[31].imag() );

		auto       fw_cpu = multi::array<complex, 1>(extensions(in_cpu));
		auto       fw_gpu = multi::thrust::cuda::array<complex, 1>(extensions(in_gpu));

		auto p_cpu = multi::fftw::plan::forward(c, in_cpu.base(), in_cpu.layout(), fw_cpu.base(), fw_cpu.layout());
		auto p_gpu = multi::cufft::plan<1>     (c,                in_gpu.layout(),                fw_gpu.layout());

		for(auto const idx : extension(in_cpu)) {
			std::cout << "B: " << idx << ": " << in_cpu[idx] << ", " << in_gpu[idx] << std::endl;
		}

		BOOST_TEST( complex(in_gpu[31]).real() == in_cpu[31].real() );
		BOOST_TEST( complex(in_gpu[31]).imag() == in_cpu[31].imag() );

		p_cpu.execute        (in_cpu.base(), fw_cpu.base());
		p_gpu.execute_forward(in_gpu.base(), fw_gpu.base());
	
		BOOST_TEST( fw_cpu[31].real() != 0.0 );
		BOOST_TEST( fw_cpu[31].imag() != 0.0 );

		for(auto const idx : extension(in_cpu)) {
			std::cout << "C: " << idx << ": " << in_cpu[idx] << ", " << in_gpu[idx] << std::endl;
		}

		BOOST_TEST( complex(in_gpu[31]).real() == in_cpu[31].real() );
		BOOST_TEST( complex(in_gpu[31]).imag() == in_cpu[31].imag() );

		for(auto const idx : extension(in_cpu)) {
			std::cout << idx << ": " << fw_cpu[idx] << ", " << fw_gpu[idx] << std::endl;
		}

		BOOST_TEST( complex(fw_gpu[31]).real() == fw_cpu[31].real() );
		BOOST_TEST( complex(fw_gpu[31]).imag() == fw_cpu[31].imag() );
	}
}

BOOST_AUTO_TEST_CASE(cufft_2D_combinations, *boost::unit_test::tolerance(0.0001)){

	using complex = thrust::complex<double>;  // this can't be std::complex<double> in the gpu

	auto const in_cpu = std::invoke([]{
		multi::array<complex, 2> ret({10, 20});
		std::default_random_engine generator;
		std::uniform_real_distribution<double> distribution(-1.0, 1.0);

		std::generate(
			reinterpret_cast<double*>(ret.data_elements()),
			reinterpret_cast<double*>(ret.data_elements() + ret.num_elements()), [&]{return distribution(generator);}
		);
		return ret;
	});

	for(auto c : std::vector<std::array<bool, 2>>{
		{true , true },
		{true , false},
		{false, true }//,
	//  {false, false}
	}){
		auto       fw_cpu = multi::array<complex, 2>(extensions(in_cpu));
		multi::fftw::dft(c, in_cpu, fw_cpu, multi::fftw::forward);

		auto const in_gpu = multi::thrust::cuda::array<complex, 2>{in_cpu};
		auto       fw_gpu = multi::thrust::cuda::array<complex, 2>(extensions(in_gpu));

		BOOST_TEST( fw_cpu[2][1].real() != 0.0 );
		BOOST_TEST( fw_cpu[2][1].imag() != 0.0 );

		multi::cufft::plan<2>(c, in_gpu.layout(), fw_gpu.layout())
			.execute(in_gpu.base(), fw_gpu.base(), multi::cufft::forward);

		BOOST_TEST( (complex(fw_gpu[2][1]) - fw_cpu[2][1]).real() == 0.0 );
		BOOST_TEST( (complex(fw_gpu[2][1]) - fw_cpu[2][1]).imag() == 0.0 );
	}
}

BOOST_AUTO_TEST_CASE(cufft_2D_combinations_inplace, *boost::unit_test::tolerance(0.0001)){

	using complex = thrust::complex<double>;  // this can't be std::complex<double> in the gpu

	auto const in_cpu = std::invoke([]{
		multi::array<complex, 2> ret({10, 20});
		std::default_random_engine generator;
		std::uniform_real_distribution<double> distribution(-1.0, 1.0);

		std::generate(
			reinterpret_cast<double*>(ret.data_elements()),
			reinterpret_cast<double*>(ret.data_elements() + ret.num_elements()), [&]{return distribution(generator);}
		);
		return ret;
	});

	for(auto c : std::vector<std::array<bool, 2>>{
		{true , true },
		{true , false},
		{false, true }//,
	//  {false, false}
	}){
		auto       fw_cpu = in_cpu;
		auto const in_gpu = multi::thrust::cuda::array<complex, 2>{in_cpu};

		multi::fftw::dft(c, fw_cpu, multi::fftw::forward);

		auto       fw_gpu = in_gpu;

		BOOST_TEST( fw_cpu[2][1].real() != 0.0 );
		BOOST_TEST( fw_cpu[2][1].imag() != 0.0 );

		multi::cufft::plan<2>(c, fw_gpu.layout(), fw_gpu.layout())
			.execute(fw_gpu.base(), fw_gpu.base(), multi::cufft::forward);

		BOOST_TEST( (complex(fw_gpu[2][1]) - fw_cpu[2][1]).real() == 0.0 );
		BOOST_TEST( (complex(fw_gpu[2][1]) - fw_cpu[2][1]).imag() == 0.0 );
	}
}

BOOST_AUTO_TEST_CASE(cufft_3D, *boost::unit_test::tolerance(0.0001)){

	using complex = thrust::complex<double>;  // this can't be std::complex<double> in the gpu

	auto const in_cpu = std::invoke([]{
		multi::array<complex, 3> ret({10, 20, 30});
		std::default_random_engine generator;
		std::uniform_real_distribution<double> distribution(-1.0, 1.0);

		std::generate(
			reinterpret_cast<double*>(ret.data_elements()), 
			reinterpret_cast<double*>(ret.data_elements() + ret.num_elements()), [&]{return distribution(generator);}
		);
		return ret;
	});

	for(auto c : std::vector<std::array<bool, 3>>{
		{true , true , true },
		{true , true , false},
		{true , false, true },
		{true , false, false},
		{false, true , true },
		{false, true , false},
		{false, false, true }//,
	//  {false, false, false}
	}){
		auto       fw_cpu = multi::array<complex, 3>(extensions(in_cpu));
		auto const in_gpu = multi::thrust::cuda::array<complex, 3>{in_cpu};

		multi::fftw::dft(c, in_cpu, fw_cpu, multi::fftw::forward);
		auto       fw_gpu = multi::thrust::cuda::array<complex, 3>(extensions(in_gpu));

		multi::cufft::dft(c, in_gpu, fw_gpu, multi::cufft::forward);

		BOOST_TEST( fw_cpu[3][2][1].real() != 0.0 );
		BOOST_TEST( fw_cpu[3][2][1].imag() != 0.0 );

		BOOST_TEST( (complex(fw_gpu[3][2][1]) - fw_cpu[3][2][1]).real() == 0.0 );
		BOOST_TEST( (complex(fw_gpu[3][2][1]) - fw_cpu[3][2][1]).imag() == 0.0 );
	}
}

BOOST_AUTO_TEST_CASE(cufft_3D_inplace, *boost::unit_test::tolerance(0.0001)){

	using complex = thrust::complex<double>;  // this can't be std::complex<double> in the gpu

	auto const in_cpu = std::invoke([]{
		multi::array<complex, 3> ret({10, 20, 30});
		std::default_random_engine generator;
		std::uniform_real_distribution<double> distribution(-1.0, 1.0);

		std::generate(
			reinterpret_cast<double*>(ret.data_elements()), 
			reinterpret_cast<double*>(ret.data_elements() + ret.num_elements()), [&]{return distribution(generator);}
		);
		return ret;
	});

	for(auto c : std::vector<std::array<bool, 3>>{
		{true , true , true },
		{true , true , false},
		{true , false, true },
		{true , false, false},
		{false, true , true },
		{false, true , false},
		{false, false, true }//,
	//  {false, false, false}
	}){
		auto       fw_cpu = in_cpu;
		auto const in_gpu = multi::thrust::cuda::array<complex, 3>{in_cpu};

		multi::fftw::dft(c, fw_cpu, multi::fftw::forward);
		auto       fw_gpu = in_gpu;

		multi::cufft::plan<3>(c, fw_gpu.layout(), fw_gpu.layout())
			.execute(fw_gpu.base(), fw_gpu.base(), multi::cufft::forward);

		BOOST_TEST( fw_cpu[3][2][1].real() != 0.0 );
		BOOST_TEST( fw_cpu[3][2][1].imag() != 0.0 );

		std::cerr << "case " << c[0] << " " << c[1] << " " << c[2] << std::endl;

		BOOST_TEST( (complex(fw_gpu[3][2][1]) - fw_cpu[3][2][1]).real() == 0.0 );
		BOOST_TEST( (complex(fw_gpu[3][2][1]) - fw_cpu[3][2][1]).imag() == 0.0 );
	}
}

BOOST_AUTO_TEST_CASE(cufft_4D, *boost::unit_test::tolerance(0.0001)){

	using complex = thrust::complex<double>;  // this can't be std::complex<double> in the gpu

	auto const in_cpu = std::invoke([]{
		multi::array<complex, 4> ret({10, 20, 30, 40});
		std::default_random_engine generator;
		std::uniform_real_distribution<double> distribution(-1.0, 1.0);

		std::generate(
			reinterpret_cast<double*>(ret.data_elements()), 
			reinterpret_cast<double*>(ret.data_elements() + ret.num_elements()), [&]{return distribution(generator);}
		);
		return ret;
	});

	for(auto c : std::vector<std::array<bool, 4>>{
		// {true , true , true , true },
		{true , true , true , false},
		{true , true , false, true },
		{true , true , false, false},
		{true , false, true , true },
		{true , false, true , false},
		{true , false, false, true },
		{true , false, false, false},
		{false, true , true , true },
		{false, true , true , false},
		{false, true , false, true },
		{false, true , false, false},
		{false, false, true , true },
		{false, false, true , false},
		{false, false, false, true }//,
	//  {false, false, false, false}
	}){
		auto       fw_cpu = multi::array<complex, 4>(extensions(in_cpu));
		multi::fftw::dft(c, in_cpu, fw_cpu, multi::fftw::forward);

		auto const in_gpu = multi::thrust::cuda::array<complex, 4>{in_cpu};
		auto       fw_gpu = multi::thrust::cuda::array<complex, 4>(extensions(in_gpu));

		BOOST_TEST( fw_cpu[4][3][2][1].real() != 0.0 );
		BOOST_TEST( fw_cpu[4][3][2][1].imag() != 0.0 );

		multi::cufft::plan<4>(c, in_gpu.layout(), fw_gpu.layout())
			.execute(in_gpu.base(), fw_gpu.base(), multi::cufft::forward);

		BOOST_TEST( (complex(fw_gpu[4][3][2][1]) - fw_cpu[4][3][2][1]).real() == 0.0 );
		BOOST_TEST( (complex(fw_gpu[4][3][2][1]) - fw_cpu[4][3][2][1]).imag() == 0.0 );
	}
}

BOOST_AUTO_TEST_CASE(cufft_3D_timing, *boost::unit_test::tolerance(0.0001)){

	auto x = multi::extensions_t<3>{300, 300, 300};
	{
		auto const in_cpu = multi::array<complex, 3>(x, 10.0);
		BOOST_ASSERT( in_cpu.num_elements()*sizeof(complex) < 2e9 );
		auto       fw_cpu = multi::array<complex, 3>(extensions(in_cpu), 99.0);
		{
		//  boost::timer::auto_cpu_timer t;  // 1.041691s wall, 1.030000s user + 0.000000s system = 1.030000s CPU (98.9%)
			multi::fftw::dft_forward({true, true}, in_cpu, fw_cpu);
			BOOST_TEST( fw_cpu[8][9][10] != 99.0 );
		}

		auto const in_gpu = multi::thrust::cuda::array<complex, 3>{in_cpu};  // (x, 10.0);
		cudaDeviceSynchronize();
		{
			auto       fw_gpu = multi::thrust::cuda::array<complex, 3>(extensions(in_gpu), 99.0);
			cudaDeviceSynchronize();
		//  boost::timer::auto_cpu_timer t; //  0.208237s wall, 0.200000s user + 0.010000s system = 0.210000s CPU (100.8%)
			boost::multi::cufft::dft({true, true}, in_gpu, fw_gpu, multi::cufft::forward);
			cudaDeviceSynchronize();
			BOOST_TEST( (static_cast<complex>(fw_gpu[8][9][10]) - fw_cpu[8][9][10]).real() == 0.0 );
			BOOST_TEST( (static_cast<complex>(fw_gpu[8][9][10]) - fw_cpu[8][9][10]).imag() == 0.0 );
		}
		{
		//  boost::timer::auto_cpu_timer t; //  0.208237s wall, 0.200000s user + 0.010000s system = 0.210000s CPU (100.8%)
			auto const fw_gpu2 = boost::multi::cufft::dft({true, true}, in_gpu, multi::cufft::forward);
			cudaDeviceSynchronize();
			BOOST_TEST( (static_cast<complex>(fw_gpu2[8][9][10]) - fw_cpu[8][9][10]).real() == 0.0 );
			BOOST_TEST( (static_cast<complex>(fw_gpu2[8][9][10]) - fw_cpu[8][9][10]).imag() == 0.0 );
		}
	}

#if 1
	{
		multi::thrust::cuda::universal_array<complex, 3> const in_gpu(x, 10.); 
		multi::thrust::cuda::universal_array<complex, 3> fw_gpu(extensions(in_gpu), 99.);

		// multi::cuda::managed::array<complex, 3> const in_gpu(x, 10.); 
		// multi::cuda::managed::array<complex, 3> fw_gpu(extensions(in_gpu), 99.);
		{
		//  boost::timer::auto_cpu_timer t; //  0.208237s wall, 0.200000s user + 0.010000s system = 0.210000s CPU (100.8%)
			multi::cufft::dft({true, true}, in_gpu, fw_gpu, multi::cufft::forward);
		//  BOOST_TEST( fw_gpu[8][9][10].operator complex() != 99. );
		}
		{
		//  boost::timer::auto_cpu_timer t; //  0.208237s wall, 0.200000s user + 0.010000s system = 0.210000s CPU (100.8%)
			multi::cufft::dft({true, true}, in_gpu, fw_gpu, multi::cufft::forward);
		//  BOOST_TEST( fw_gpu[8][9][10].operator complex() != 99. );
		}
	}
#endif
}

#if 0

BOOST_AUTO_TEST_CASE(cufft_combinations, *utf::tolerance(0.00001)){

	auto const in = []{
		multi::array<complex, 4> ret({32, 90, 98, 96});
		std::generate(ret.data_elements(), ret.data_elements() + ret.num_elements(), 
			[](){return complex{std::rand()*1./RAND_MAX, std::rand()*1./RAND_MAX};}
		);
		return ret;
	}();
	std::clog<<"memory size "<< in.num_elements()*sizeof(complex)/1e6 <<" MB\n";

	multi::thrust::cuda::universal_array<complex, 4> const in_gpu = in;
	multi::thrust::cuda::universal_array<complex, 4> const in_mng = in;

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
		//  BOOST_TEST( abs( static_cast<multi::complex<double>>(in_rw[5][4][3][1]) - multi::complex<double>(out[5][4][3][1]) ) == 0. );            
		}();
		{
			multi::array<complex, 4> in_rw2 = in;
			[&, _ = watch{"cpu_mov "}]{
				multi::array<complex, 4> const out_mov = multi::fftw::dft_forward(c, std::move(in_rw2));
			//  what(out_mov);
			//  BOOST_TEST( abs( static_cast<multi::complex<double>>(out_mov[5][4][3][1]) - multi::complex<double>(out[5][4][3][1]) ) == 0. );          
				BOOST_REQUIRE( is_empty(in_rw2) );
				BOOST_REQUIRE( extensions(out_mov) == extensions(in) );
			}();
		}


		[&, _ = watch{"cpu_new "}]{
			auto const out_cpy = multi::fftw::dft_forward(c, in);
			BOOST_TEST( abs( static_cast<std::complex<double>>(out_cpy[5][4][3][1]) - std::complex<double>(out[5][4][3][1]) ) == 0. );
		}();
		multi::thrust::cuda::array<complex, 4> out_gpu(extensions(in_gpu));
		[&, _ = watch{"gpu_opl "}]{
			multi::cufft::dft(c, in_gpu   , out_gpu, multi::cufft::forward);
			BOOST_TEST( abs( static_cast<complex>(out_gpu[5][4][3][1]) - out[5][4][3][1] ) == 0. );
		}();
		{
			multi::thrust::cuda::array<complex, 4> in_rw_gpu = in_gpu;
			[&, _ = watch{"gpu_ipl "}]{
				multi::cufft::dft(c, in_rw_gpu, multi::cufft::forward);
				BOOST_TEST( abs( static_cast<complex>(in_rw_gpu[5][4][3][1]) - out[5][4][3][1] ) == 0. );
			}();
		}
		{
			multi::thrust::cuda::array<complex, 4> in_rw_gpu = in_gpu;
			[&, _ = watch{"gpu_mov "}]{
				multi::thrust::cuda::array<complex, 4> const out_mov = multi::cufft::dft_forward(c, std::move(in_rw_gpu));
			//  BOOST_REQUIRE( in_rw_gpu.empty() );
			//  BOOST_TEST( abs( static_cast<complex>(out_mov[5][4][3][1]) - out[5][4][3][1] ) == 0. );
			}();
		}
		{
			multi::thrust::cuda::array<complex, 4> in_rw_gpu = in_gpu;
			[&, _ = watch{"gpu_mov "}]{
				multi::thrust::cuda::array<complex, 4> out_mov = std::move(in_rw_gpu);
				multi::cufft::dft(c, out_mov, multi::cufft::forward);
			//  BOOST_REQUIRE( in_rw_gpu.empty() );
			//  BOOST_TEST( abs( static_cast<complex>(out_mov[5][4][3][1]) - out[5][4][3][1] ) == 0. );
			}();
		}
		cudaDeviceSynchronize();
		[&, _ = watch{"gpu_new "}]{
			multi::thrust::cuda::array<complex, 4> const out_cpy = multi::cufft::dft(c, in_gpu, multi::cufft::forward);
		}();
		multi::thrust::cuda::universal_array<complex, 4> out_mng(extensions(in_mng));
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
	// std::clog<<"cache size "
	//  << multi::cufft::plan::cache<1>().size() <<' '
	//  << multi::cufft::plan::cache<2>().size() <<' '
	//  << multi::cufft::plan::cache<3>().size() <<' '
	//  << multi::cufft::plan::cache<4>().size() <<' '
	// <<std::endl;
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

	multi::thrust::cuda::array<complex, 4> const in = in_cpu;
	multi::thrust::cuda::array<complex, 4>       out(extensions(in));

#if 0
	multi::cufft::many_dft(begin(unrotated(in)), end(unrotated(in)), begin(unrotated(out)), +1);

	multi::array<complex, 4> out_cpu(extensions(in));
	multi::fft::many_dft(begin(unrotated(in_cpu)), end(unrotated(in_cpu)), begin(unrotated(out_cpu)), +1);

	BOOST_TEST( imag( static_cast<complex>(out[5][4][3][2]) - out_cpu[5][4][3][2]) == 0. );
#endif
}

#if 0
BOOST_AUTO_TEST_CASE(cufft_4D, *utf::tolerance(0.00001) ){
	auto const in = []{
		multi::array<complex, 3> ret({10, 10, 10});
		std::generate(ret.data_elements(), ret.data_elements() + ret.num_elements(), 
			[](){return complex{std::rand()*1./RAND_MAX, std::rand()*1./RAND_MAX};}
		);
		return ret;
	}();

	multi::array<complex, 3> out(extensions(in));
//  multi::fftw::dft({true, false, true}, in, out, multi::fftw::forward);
	multi::fftw::many_dft(begin(in.rotated()), end(in.rotated()), begin(out.rotated()), multi::fftw::forward);

	multi::thrust::cuda::array<complex, 3> in_gpu = in;
	multi::thrust::cuda::array<complex, 3> out_gpu(extensions(in));

//  multi::cufft::dft({true, false, true}, in_gpu, out_gpu, multi::fft::forward);//multi::cufft::forward);  
	// multi::cufft::many_dft(begin(in_gpu.rotated()), end(in_gpu.rotated()), begin( out_gpu.rotated() ), multi::fftw::forward);
	// BOOST_TEST( ( static_cast<complex>(out_gpu[5][4][3]) - out[5][4][3]).imag() == 0. );
}
#endif

#endif
