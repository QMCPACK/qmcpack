// Copyright 2020-2025 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/multi/adaptors/fftw.hpp>
#include <boost/multi/array.hpp>

#include <boost/core/lightweight_test.hpp>

#include <numeric>
#include <random>
#include <stdexcept>

#if(!(defined(__HIP_PLATFORM_AMD__) || defined(__HIP_PLATFORM_NVIDIA__))) && (!defined(__HIPCC__))
#include <boost/multi/adaptors/cufft.hpp>
#else
#include <boost/multi/adaptors/hipfft.hpp>
#endif

#include <boost/multi/adaptors/fft.hpp>
#include <boost/multi/adaptors/thrust.hpp>

#include <thrust/complex.h>
#include <thrust/transform_reduce.h>

#include <cassert>
#include <chrono>
#include <iostream>
#include <string>
#include <utility>

namespace multi = boost::multi;
using complex   = thrust::complex<double>;

template<>
constexpr bool multi::force_element_trivial_default_construction<thrust::complex<double>> = true;

template<class T>
__attribute__((always_inline)) inline void DoNotOptimize(T const& value) {  // NOLINT(readability-identifier-naming) consistency with Google benchmark
	asm volatile("" : "+m"(const_cast<T&>(value)));                         // NOLINT(hicpp-no-assembler,cppcoreguidelines-pro-type-const-cast) hack
}

class watch : std::chrono::high_resolution_clock {
	std::string label_;
	time_point  start_;

 public:
	explicit watch(char const* label) : label_{label} {
		cudaDeviceSynchronize() == cudaSuccess ? void() : assert(0);  // NOLINT(misc-include-cleaner) the header is included conditionally
		start_ = now();
	}

	watch(watch const&) = delete;
	watch(watch&&)      = delete;

	auto operator=(watch const&) -> watch& = delete;
	auto operator=(watch&&) -> watch&      = delete;

	watch() : watch("") {}
	~watch() {
		cudaDeviceSynchronize() == cudaSuccess ? void() : assert(0);
		auto const count = std::chrono::duration<double>(now() - start_).count();
		std::cerr << label_ << ": " << count << " sec\n";
	}
};

using complex = thrust::complex<double>;  // this can't be std::complex<double> in the gpu

struct norm_t {
	__host__ __device__ auto operator()(complex const& x) const {
		return thrust::norm(x);
	}
};

auto main() -> int try {
	complex const I{0.0, 1.0};  // NOLINT(readability-identifier-length)

	// BOOST_AUTO_TEST_CASE(cufft_2D, *boost::unit_test::tolerance(0.0001))
	{
		auto const in_cpu = multi::array<complex, 2>{
			{ 1.0 + 2.0 * I,  9.0 - 1.0 * I,  2.0 + 4.0 * I},
			{ 3.0 + 3.0 * I,  7.0 - 4.0 * I,  1.0 + 9.0 * I},
			{ 4.0 + 1.0 * I,  5.0 + 3.0 * I,  2.0 + 4.0 * I},
			{ 3.0 - 1.0 * I,  8.0 + 7.0 * I,  2.0 + 1.0 * I},
			{31.0 - 1.0 * I, 18.0 + 7.0 * I, 2.0 + 10.0 * I}
		};

		{
			auto fw_cpu = multi::array<complex, 2>(extents(in_cpu));
			multi::fftw::dft_forward({true, true}, in_cpu, fw_cpu);

			auto const in_gpu = multi::thrust::cuda::array<complex, 2>{in_cpu};
			auto       fw_gpu = multi::thrust::cuda::array<complex, 2>(extents(in_gpu));

			BOOST_TEST( fw_cpu[3][2].real() != 0.0 );
			BOOST_TEST( fw_cpu[3][2].imag() != 0.0 );

			multi::cufft::plan<2>({true, true}, in_gpu.layout(), fw_gpu.layout())
				.execute(in_gpu.base(), fw_gpu.base(), multi::cufft::forward);

			BOOST_TEST( std::abs((complex(fw_gpu[3][2]) - fw_cpu[3][2]).real()) < 1.0e-8 );
			BOOST_TEST( std::abs((complex(fw_gpu[3][2]) - fw_cpu[3][2]).imag()) < 1.0e-8 );

			// TODO(correaa) test funcional interface for GPU
			// auto const& dft = multi::fft::DFT({true, true}, in_cpu, multi::fft::forward);

			// BOOST_TEST( dft.extents() == in_cpu.extents() );
			// BOOST_TEST( (*dft.begin()).size() == (*in_cpu.begin()).size() );
			// BOOST_TEST( (*dft.begin()).extents() == (*in_cpu.begin()).extents() );

			// multi::array<complex, 2> const fw_cpu_out = multi::fft::DFT({true, true}, in_cpu, multi::fft::forward);
		}
	}
	{
		auto const in_cpu = std::invoke([] {
			multi::array<complex, 4> ret({20, 20, 20, 20});
			auto const [is, js, ks, ls] = ret.extents();
			for(auto i : is)
				for(auto j : js)
					for(auto k : ks)
						for(auto l : ls) {
							ret[i][j][k][l] = complex{
								static_cast<double>(i + j + k + l),
								static_cast<double>(i - j + k - l),
							};
						}
			return ret;
		});

		auto const in_gpu = multi::thrust::cuda::array<complex, 4>{in_cpu};

		auto const nrm = thrust::transform_reduce(
			in_gpu.elements().begin(), in_gpu.elements().end(),
			norm_t{}, 0.0, thrust::plus<>{}
		);

		auto fw_gpu = multi::thrust::cuda::array<complex, 4>(in_gpu.extents());
		fw_gpu      = in_gpu;
		// multi::cufft::plan<4>({true, true, true, true}, in_gpu.layout(), fw_gpu.layout())
		//  .execute(in_gpu.base(), fw_gpu.base(), multi::cufft::forward);

		// cudaDeviceSynchronize() == cudaSuccess ? void() : throw std::runtime_error{"cuda error"};

		// auto const nrm_fwd = thrust::transform_reduce(
		//  fw_gpu.elements().begin(), fw_gpu.elements().end(),
		//  norm_t{}, 0.0, thrust::plus<>{}
		// );
		// std::cout << "norm: " << nrm*20.0*20.0 << ", norm forward: " << nrm_fwd << '\n';
		// BOOST_TEST( nrm_fwd == nrm*20.0*20.0 );
	}
	{
		auto const in_cpu = multi::array<complex, 2>{
			{ 1.0 + 2.0 * I,  9.0 - 1.0 * I,  2.0 + 4.0 * I},
			{ 3.0 + 3.0 * I,  7.0 - 4.0 * I,  1.0 + 9.0 * I},
			{ 4.0 + 1.0 * I,  5.0 + 3.0 * I,  2.0 + 4.0 * I},
			{ 3.0 - 1.0 * I,  8.0 + 7.0 * I,  2.0 + 1.0 * I},
			{31.0 - 1.0 * I, 18.0 + 7.0 * I, 2.0 + 10.0 * I}
		};

		auto fw_cpu = multi::array<complex, 2>(extents(in_cpu));
		multi::fftw::dft({false, true}, in_cpu, fw_cpu, multi::fftw::forward);

		auto const in_gpu = multi::thrust::cuda::array<complex, 2>{in_cpu};
		auto       fw_gpu = multi::thrust::cuda::array<complex, 2>(extents(in_gpu));

		BOOST_TEST( fw_cpu[3][2].real() != 0.0 );
		BOOST_TEST( fw_cpu[3][2].imag() != 0.0 );

		multi::cufft::plan<2>({false, true}, in_gpu.layout(), fw_gpu.layout())
			.execute(in_gpu.base(), fw_gpu.base(), multi::cufft::forward);

		BOOST_TEST( thrust::abs(complex(fw_gpu[3][2]) - fw_cpu[3][2]) <  1e-12 );
	}
	{
		auto const in_cpu = multi::array<complex, 2>{
			{ 1.0 + 2.0 * I,  9.0 - 1.0 * I,  2.0 + 4.0 * I},
			{ 3.0 + 3.0 * I,  7.0 - 4.0 * I,  1.0 + 9.0 * I},
			{ 4.0 + 1.0 * I,  5.0 + 3.0 * I,  2.0 + 4.0 * I},
			{ 3.0 - 1.0 * I,  8.0 + 7.0 * I,  2.0 + 1.0 * I},
			{31.0 - 1.0 * I, 18.0 + 7.0 * I, 2.0 + 10.0 * I}
		};

		auto fw_cpu = multi::array<complex, 2>(extents(in_cpu));
		multi::fftw::dft({false, true}, in_cpu, fw_cpu, multi::fftw::forward);

		auto const in_gpu = multi::thrust::cuda::array<complex, 2>{in_cpu};
		auto       fw_gpu = multi::thrust::cuda::array<complex, 2>(extents(in_gpu));

		BOOST_TEST( fw_cpu[3][2].real() != 0.0 );
		BOOST_TEST( fw_cpu[3][2].imag() != 0.0 );

		for(int i = 0; i != in_gpu.size(); ++i) {
			multi::cufft::plan<1>({true}, in_gpu[i].layout(), fw_gpu[i].layout())
				.execute(in_gpu[i].base(), fw_gpu[i].base(), multi::cufft::forward);
		}

		BOOST_TEST( thrust::abs(complex(fw_gpu[3][2]) - fw_cpu[3][2]) < 1.0e-12 );
	}
	{
		auto const in_cpu = multi::array<complex, 2>{
			{ 1.0 + 2.0 * I,  9.0 - 1.0 * I,  2.0 + 4.0 * I},
			{ 3.0 + 3.0 * I,  7.0 - 4.0 * I,  1.0 + 9.0 * I},
			{ 4.0 + 1.0 * I,  5.0 + 3.0 * I,  2.0 + 4.0 * I},
			{ 3.0 - 1.0 * I,  8.0 + 7.0 * I,  2.0 + 1.0 * I},
			{31.0 - 1.0 * I, 18.0 + 7.0 * I, 2.0 + 10.0 * I}
		};

		auto fw_cpu = multi::array<complex, 2>(extents(in_cpu));
		multi::fftw::dft({false, true}, in_cpu, fw_cpu, multi::fftw::forward);

		auto const in_gpu  = multi::thrust::cuda::array<complex, 2>{in_cpu};
		auto       fw_gpu  = multi::thrust::cuda::array<complex, 2>(extents(in_gpu));
		auto       fw_gpu2 = multi::thrust::cuda::array<complex, 2>(extents(in_gpu));
		auto       fw_gpu3 = multi::thrust::cuda::array<complex, 2>(extents(in_gpu));

		BOOST_TEST( fw_cpu[3][2].real() != 0.0 );
		BOOST_TEST( fw_cpu[3][2].imag() != 0.0 );

		for(int i = 0; i != in_gpu.size(); ++i) {
			multi::cufft::plan<1>({true}, in_gpu[i].layout(), fw_gpu[i].layout())
				.execute(in_gpu[i].base(), fw_gpu[i].base(), multi::cufft::forward);
		}

		multi::cufft::plan<2>({false, true}, in_gpu.layout(), fw_gpu2.layout())
			.execute(in_gpu.base(), fw_gpu2.base(), multi::cufft::forward);

		BOOST_TEST( abs(complex(fw_gpu[3][2]) - fw_cpu[3][2]) < 1e-10 );
		BOOST_TEST( abs(complex(fw_gpu[3][2]) - complex(fw_gpu2[3][2])) < 1e-10 );
	}

	{
		auto const in_cpu = multi::array<complex, 2>{
			{ 1.0 + 2.0 * I,  9.0 - 1.0 * I,  2.0 + 4.0 * I},
			{ 3.0 + 3.0 * I,  7.0 - 4.0 * I,  1.0 + 9.0 * I},
			{ 4.0 + 1.0 * I,  5.0 + 3.0 * I,  2.0 + 4.0 * I},
			{ 3.0 - 1.0 * I,  8.0 + 7.0 * I,  2.0 + 1.0 * I},
			{31.0 - 1.0 * I, 18.0 + 7.0 * I, 2.0 + 10.0 * I}
		};
		auto fw_cpu = multi::array<complex, 2>(extents(in_cpu));
		multi::fftw::dft({false, true}, in_cpu, fw_cpu, multi::fftw::forward);

		auto const in_gpu = multi::thrust::cuda::array<complex, 2>{in_cpu};
		auto const fw_gpu = multi::cufft::dft({false, true}, in_gpu, multi::cufft::forward);

		BOOST_TEST( abs(fw_cpu[3][2]) != 0.0 );

		BOOST_TEST( abs(complex(fw_gpu[3][2]) - fw_cpu[3][2]) < 1e-10 );
	}
	{
		auto const in_cpu = multi::array<complex, 2>{
			{ 1.0 + 2.0 * I,  9.0 - 1.0 * I,  2.0 + 4.0 * I},
			{ 3.0 + 3.0 * I,  7.0 - 4.0 * I,  1.0 + 9.0 * I},
			{ 4.0 + 1.0 * I,  5.0 + 3.0 * I,  2.0 + 4.0 * I},
			{ 3.0 - 1.0 * I,  8.0 + 7.0 * I,  2.0 + 1.0 * I},
			{31.0 - 1.0 * I, 18.0 + 7.0 * I, 2.0 + 10.0 * I}
		};
		auto fw_cpu = multi::array<complex, 2>(extents(in_cpu));
		multi::fftw::dft({true, false}, in_cpu, fw_cpu, multi::fftw::forward);

		auto const in_gpu = multi::thrust::cuda::array<complex, 2>{in_cpu};
		auto const fw_gpu = multi::cufft::dft({true, false}, in_gpu, multi::cufft::forward);

		BOOST_TEST( fw_cpu.extents() == in_cpu.extents() );
		BOOST_TEST( abs(fw_cpu[3][2]) != 0.0 );

		BOOST_TEST( fw_gpu.extents() == in_gpu.extents() );
		BOOST_TEST( abs(complex(fw_gpu[3][2]) - fw_cpu[3][2]) < 1e-10 );
		BOOST_TEST( abs(complex(fw_gpu[2][1]) - fw_cpu[2][1]) < 1e-10 );
	}

	// BOOST_AUTO_TEST_CASE(cufft_1D_combinations, *boost::unit_test::tolerance(0.0001))
	{
		using complex = thrust::complex<double>;  // this can't be std::complex<double> in the gpu

		auto const in_cpu = std::invoke([] {
			multi::array<complex, 1>               ret({128}, complex{});
			std::default_random_engine             generator;
			std::uniform_real_distribution<double> distribution(1.0, 88.0);

			std::generate(
				reinterpret_cast<double*>(ret.data_elements()),
				reinterpret_cast<double*>(ret.data_elements() + ret.num_elements()), [&] { return distribution(generator); }
			);
			return ret;
		});

		for(auto c : std::vector<std::array<bool, 1>>{
				{true}  //,
						// {false},
			}) {
			auto const in_gpu = multi::thrust::cuda::array<complex, 1>{in_cpu};

			BOOST_TEST( complex(in_gpu[31]).real() == in_cpu[31].real() );
			BOOST_TEST( complex(in_gpu[31]).imag() == in_cpu[31].imag() );

			auto fw_cpu = multi::array<complex, 1>(extents(in_cpu));
			auto fw_gpu = multi::thrust::cuda::array<complex, 1>(extents(in_gpu));

			auto p_cpu = multi::fftw::plan::forward(c, in_cpu.base(), in_cpu.layout(), fw_cpu.base(), fw_cpu.layout());
			auto p_gpu = multi::cufft::plan<1>(c, in_gpu.layout(), fw_gpu.layout());

			BOOST_TEST( abs(complex(in_gpu[31]) -  in_cpu[31]) < 1e-10 );

			p_cpu.execute(in_cpu.base(), fw_cpu.base());
			p_gpu.execute_forward(in_gpu.base(), fw_gpu.base());

			BOOST_TEST( abs(fw_cpu[31]) != 0.0 );

			BOOST_TEST( abs( complex(in_gpu[31]) - in_cpu[31]) < 1e-10 );
			BOOST_TEST( abs( complex(fw_gpu[31]) - fw_cpu[31]) < 1e-10 );
		}
	}

	// BOOST_AUTO_TEST_CASE(cufft_2D_combinations, *boost::unit_test::tolerance(0.0001))
	{

		using complex = thrust::complex<double>;  // this can't be std::complex<double> in the gpu

		auto const in_cpu = std::invoke([] {
			multi::array<complex, 2>               ret({10, 20});
			std::default_random_engine             generator;
			std::uniform_real_distribution<double> distribution(-1.0, 1.0);

			std::generate(
				reinterpret_cast<double*>(ret.data_elements()),
				reinterpret_cast<double*>(ret.data_elements() + ret.num_elements()), [&] { return distribution(generator); }
			);
			return ret;
		});

		for(auto c : std::vector<std::array<bool, 2>>{
				{ true,  true},
				{ true, false},
				{false,  true}, //  {false, false}
        }) {
			auto fw_cpu = multi::array<complex, 2>(extents(in_cpu));
			multi::fftw::dft(c, in_cpu, fw_cpu, multi::fftw::forward);

			auto const in_gpu = multi::thrust::cuda::array<complex, 2>{in_cpu};
			auto       fw_gpu = multi::thrust::cuda::array<complex, 2>(extents(in_gpu));

			BOOST_TEST( abs(fw_cpu[2][1]) != 0.0 );

			multi::cufft::plan<2>(c, in_gpu.layout(), fw_gpu.layout())
				.execute(in_gpu.base(), fw_gpu.base(), multi::cufft::forward);

			BOOST_TEST( abs(complex(fw_gpu[2][1]) - fw_cpu[2][1]) < 1e-10 );
		}
	}

	// BOOST_AUTO_TEST_CASE(cufft_2D_combinations_inplace, *boost::unit_test::tolerance(0.0001))
	{

		using complex = thrust::complex<double>;  // this can't be std::complex<double> in the gpu

		auto const in_cpu = std::invoke([] {
			multi::array<complex, 2>               ret({10, 20});
			std::default_random_engine             generator;
			std::uniform_real_distribution<double> distribution(-1.0, 1.0);

			std::generate(
				reinterpret_cast<double*>(ret.data_elements()),
				reinterpret_cast<double*>(ret.data_elements() + ret.num_elements()), [&] { return distribution(generator); }
			);
			return ret;
		});

		for(auto c : std::vector<std::array<bool, 2>>{
				{ true,  true},
				{ true, false},
				{false,  true}  //,
							   //  {false, false}
        }) {
			auto       fw_cpu = in_cpu;
			auto const in_gpu = multi::thrust::cuda::array<complex, 2>{in_cpu};

			multi::fftw::dft(c, fw_cpu, multi::fftw::forward);

			auto fw_gpu = in_gpu;

			BOOST_TEST( abs(fw_cpu[2][1]) != 0.0 );

			multi::cufft::plan<2>(c, fw_gpu.layout(), fw_gpu.layout())
				.execute(fw_gpu.base(), fw_gpu.base(), multi::cufft::forward);

			BOOST_TEST( abs(complex(fw_gpu[2][1]) - fw_cpu[2][1]) < 1e-10 );
		}
	}

	// BOOST_AUTO_TEST_CASE(cufft_3D, *boost::unit_test::tolerance(0.0001))
	{

		using complex = thrust::complex<double>;  // this can't be std::complex<double> in the gpu

		auto const in_cpu = std::invoke([] {
			multi::array<complex, 3>               ret({10, 20, 30});
			std::default_random_engine             generator;
			std::uniform_real_distribution<double> distribution(-1.0, 1.0);

			std::generate(
				reinterpret_cast<double*>(ret.data_elements()),
				reinterpret_cast<double*>(ret.data_elements() + ret.num_elements()), [&] { return distribution(generator); }
			);
			return ret;
		});

		for(auto c : std::vector<std::array<bool, 3>>{
				{ true,  true,  true},
				{ true,  true, false},
				{ true, false,  true},
				{ true, false, false},
				{false,  true,  true},
				{false,  true, false},
				{false, false,  true}  //,
									  //  {false, false, false}
        }) {
			auto       fw_cpu = multi::array<complex, 3>(extents(in_cpu));
			auto const in_gpu = multi::thrust::cuda::array<complex, 3>{in_cpu};

			multi::fftw::dft(c, in_cpu, fw_cpu, multi::fftw::forward);
			auto fw_gpu = multi::thrust::cuda::array<complex, 3>(extents(in_gpu));

			multi::cufft::dft(c, in_gpu, fw_gpu, multi::cufft::forward);

			BOOST_TEST( abs(fw_cpu[3][2][1]) != 0.0 );

			BOOST_TEST( abs(complex(fw_gpu[3][2][1]) - fw_cpu[3][2][1]) < 1e-10 );
		}
	}

	// BOOST_AUTO_TEST_CASE(cufft_3D_inplace, *boost::unit_test::tolerance(0.0001))
	{

		using complex = thrust::complex<double>;  // this can't be std::complex<double> in the gpu

		auto const in_cpu = std::invoke([] {
			multi::array<complex, 3>               ret({10, 20, 30});
			std::default_random_engine             generator;
			std::uniform_real_distribution<double> distribution(-1.0, 1.0);

			std::generate(
				reinterpret_cast<double*>(ret.data_elements()),
				reinterpret_cast<double*>(ret.data_elements() + ret.num_elements()), [&] { return distribution(generator); }
			);
			return ret;
		});

		for(auto c : std::vector<std::array<bool, 3>>{
				{ true,  true,  true},
				{ true,  true, false},
				{ true, false,  true},
				{ true, false, false},
				{false,  true,  true},
				{false,  true, false},
				{false, false,  true}  //,
									  //  {false, false, false}
        }) {
			auto       fw_cpu = in_cpu;
			auto const in_gpu = multi::thrust::cuda::array<complex, 3>{in_cpu};

			multi::fftw::dft(c, fw_cpu, multi::fftw::forward);
			auto fw_gpu = in_gpu;

			multi::cufft::plan<3>(c, fw_gpu.layout(), fw_gpu.layout())
				.execute(fw_gpu.base(), fw_gpu.base(), multi::cufft::forward);

			BOOST_TEST( abs(fw_cpu[3][2][1]) != 0.0 );

			// std::cerr << "case " << c[0] << " " << c[1] << " " << c[2] << std::endl;
			// std::cerr << complex(fw_gpu[3][2][1]) - fw_cpu[3][2][1] << std::endl;
			// BOOST_TEST( abs(complex(fw_gpu[3][2][1]) - fw_cpu[3][2][1]) < 1e-10 );
			// TODO(correaa), these two cases are failing
			// case 1 1 1 * (-34.154,-39.0958)
			// case 1 1 0   (0,-1.77636e-15)
			// case 1 0 1 * (-12.6338,0.299744)
			// case 1 0 0 * (4.44089e-16,-4.44089e-16)
			// case 0 1 1   (20.1121,-10.8888)
			// case 0 1 0 * (0,-2.22045e-16)
			// case 0 0 1   (-0.348103,4.32914)
		}
	}

	// BOOST_AUTO_TEST_CASE(cufft_4D, *boost::unit_test::tolerance(0.0001)
	{

		using complex = thrust::complex<double>;  // this can't be std::complex<double> in the gpu

		auto const in_cpu = std::invoke([] {
			multi::array<complex, 4>               ret({10, 20, 30, 40});
			std::default_random_engine             generator;
			std::uniform_real_distribution<double> distribution(-1.0, 1.0);

			std::generate(
				reinterpret_cast<double*>(ret.data_elements()),
				reinterpret_cast<double*>(ret.data_elements() + ret.num_elements()), [&] { return distribution(generator); }
			);
			return ret;
		});

		for(auto c : std::vector<std::array<bool, 4>>{
				{true , true , true , true },
				{ true,  true,  true, false},
				{ true,  true, false,  true},
				{ true,  true, false, false},
				{ true, false,  true,  true},
				{ true, false,  true, false},
				{ true, false, false,  true},
				{ true, false, false, false},
				{false,  true,  true,  true},
				{false,  true,  true, false},
				{false,  true, false,  true},
				{false,  true, false, false},
				{false, false,  true,  true},
				{false, false,  true, false},
				{false, false, false,  true}  //,
				//  {false, false, false, false}
        }) {
			auto fw_cpu = multi::array<complex, 4>(extents(in_cpu));
			multi::fftw::dft(c, in_cpu, fw_cpu, multi::fftw::forward);

			auto const in_gpu = multi::thrust::cuda::array<complex, 4>{in_cpu};
			auto       fw_gpu = multi::thrust::cuda::array<complex, 4>(extents(in_gpu));

			BOOST_TEST( abs(fw_cpu[4][3][2][1]) != 0 );

			multi::cufft::dft(c, in_gpu, fw_gpu, multi::cufft::forward);

			std::cerr << "Case " << c[0] << " " << c[1] << " " << c[2] << " " << c[3] << ": " << complex(fw_gpu[4][3][2][1]) - fw_cpu[4][3][2][1] << '\n';

			BOOST_TEST( abs(complex(fw_gpu[4][3][2][1]) - fw_cpu[4][3][2][1]) < 1e-10 );
		}
	}

	// SO 4D intermediate answer: 
	{
		using complex = thrust::complex<double>;  // this can't be std::complex<double> in the gpu

		auto const in_cpu = std::invoke([] {
			multi::array<complex, 4>               ret({12, 128, 128, 4});
			std::default_random_engine             generator;
			std::uniform_real_distribution<double> distribution(-1.0, 1.0);

			std::generate(
				ret.elements().begin(), ret.elements().end(), [&] { return distribution(generator); }
			);
			return ret;
		});

		auto fw_cpu = multi::array<complex, 4>(extents(in_cpu));
		multi::fftw::dft({false, true, true, false}, in_cpu, fw_cpu, multi::fftw::forward);

		auto const in_gpu = multi::thrust::cuda::array<complex, 4>{in_cpu};
		auto       fw_gpu = multi::thrust::cuda::array<complex, 4>(extents(in_gpu));

		BOOST_TEST( abs(fw_cpu[4][3][2][1]) != 0.0 );

		multi::cufft::dft({false, true, true, false}, in_gpu, fw_gpu, multi::cufft::forward);

		BOOST_TEST( abs(complex(fw_gpu[4][3][2][1]) - fw_cpu[4][3][2][1]) < 1e-10 );
	}
	// small case 99
	{
		using complex = thrust::complex<double>;  // this can't be std::complex<double> in the gpu

		auto const in_cpu = std::invoke([] {
			multi::array<complex, 4>               ret({12, 128, 128, 4});
			std::default_random_engine             generator;
			std::uniform_real_distribution<double> distribution(-1.0, 1.0);

			std::generate(
				ret.elements().begin(), ret.elements().end(), [&] { return distribution(generator); }
			);
			return ret;
		});

		multi::thrust::cuda::array<complex, 4>       in({12, 128, 128, 4});
		in = in_cpu;

		multi::thrust::cuda::array<complex, 4>       ou({12, 128, 128, 4}, 0.0);

		multi::cufft::dft_forward({false, true, true, false}, in, ou);

		std::cout << "small case : " << ou[4][3][2][1] << '\n';
	}

	return boost::report_errors();
} catch(...) {
	throw;
	return 1;
}

// #if 0

// }

// BOOST_AUTO_TEST_CASE(check_thrust_complex_vs_std_complex, *boost::unit_test::tolerance(0.0001)){

//  multi::array<std   ::complex<double>, 1> const s_in = {1.0 + I*2.0, 2.0 + I*3.0, 3.0 + I*4.0};
//  multi::array<thrust::complex<double>, 1> const t_in = {1.0 + I*2.0, 2.0 + I*3.0, 3.0 + I*4.0};

//  multi::array<std   ::complex<double>, 1>       s_out(s_in.extents());
//  multi::array<thrust::complex<double>, 1>       t_out(t_in.extents());

//  multi::fftw::plan::forward({true}, s_in.base(), s_in.layout(), s_out.base(), s_out.layout()).execute(s_in.base(), s_out.base());
//  multi::fftw::plan::forward({true}, t_in.base(), t_in.layout(), t_out.base(), t_out.layout()).execute(t_in.base(), t_out.base());

//  BOOST_REQUIRE( std::equal(s_out.begin(), s_out.end(), t_out.begin()) );
// }

// BOOST_AUTO_TEST_CASE(small_1D_cpu_vs_cpu, *boost::unit_test::tolerance(0.0001)){

//  multi::array<thrust::complex<double>, 1> const cpu_in = {1.0 + I*2.0, 2.0 + I*3.0, 3.0 + I*4.0};
//  multi::thrust::cuda::array<thrust::complex<double>, 1> const gpu_in = {1.0 + I*2.0, 2.0 + I*3.0, 3.0 + I*4.0};

//  multi::array<thrust::complex<double>, 1> cpu_out(cpu_in.extents());
//  multi::thrust::cuda::array<thrust::complex<double>, 1> gpu_out(gpu_in.extents());

//  multi::fftw::plan::forward({true}, cpu_in.base(), cpu_in.layout(), cpu_out.base(), cpu_out.layout()).execute        (cpu_in.base(), cpu_out.base());
//  multi::cufft::plan<1>     ({true},                gpu_in.layout(),                 gpu_out.layout()).execute_forward(gpu_in.base(), gpu_out.base());
// }

// BOOST_AUTO_TEST_CASE(cufft_3D_timing, *boost::unit_test::tolerance(0.0001)){

//  auto x = multi::extents_t<3>{300, 300, 300};
//  {
//      auto const in_cpu = multi::array<complex, 3>(x, 10.0);
//      BOOST_ASSERT( in_cpu.num_elements()*sizeof(complex) < 2e9 );
//      auto       fw_cpu = multi::array<complex, 3>(extents(in_cpu), 99.0);
//      {
//      //  boost::timer::auto_cpu_timer t;  // 1.041691s wall, 1.030000s user + 0.000000s system = 1.030000s CPU (98.9%)
//          multi::fftw::dft_forward({true, true}, in_cpu, fw_cpu);
//          BOOST_TEST( fw_cpu[8][9][10] != 99.0 );
//      }

//      auto const in_gpu = multi::thrust::cuda::array<complex, 3>{in_cpu};  // (x, 10.0);
//      cudaDeviceSynchronize()==cudaSuccess?void():assert(0);
//      {
//          auto       fw_gpu = multi::thrust::cuda::array<complex, 3>(extents(in_gpu), 99.0);
//          cudaDeviceSynchronize()==cudaSuccess?void():assert(0);
//      //  boost::timer::auto_cpu_timer t; //  0.208237s wall, 0.200000s user + 0.010000s system = 0.210000s CPU (100.8%)
//          boost::multi::cufft::dft({true, true}, in_gpu, fw_gpu, multi::cufft::forward);
//          cudaDeviceSynchronize()==cudaSuccess?void():assert(0);
//          BOOST_TEST( (static_cast<complex>(fw_gpu[8][9][10]) - fw_cpu[8][9][10]).real() == 0.0 );
//          BOOST_TEST( (static_cast<complex>(fw_gpu[8][9][10]) - fw_cpu[8][9][10]).imag() == 0.0 );
//      }
//      {
//      //  boost::timer::auto_cpu_timer t; //  0.208237s wall, 0.200000s user + 0.010000s system = 0.210000s CPU (100.8%)
//          auto const fw_gpu2 = boost::multi::cufft::dft({true, true}, in_gpu, multi::cufft::forward);
//          cudaDeviceSynchronize()==cudaSuccess?void():assert(0);
//          BOOST_TEST( (static_cast<complex>(fw_gpu2[8][9][10]) - fw_cpu[8][9][10]).real() == 0.0 );
//          BOOST_TEST( (static_cast<complex>(fw_gpu2[8][9][10]) - fw_cpu[8][9][10]).imag() == 0.0 );
//      }
//  }

// #if 1
//  {
//      multi::thrust::cuda::universal_array<complex, 3> const in_gpu(x, 10.);
//      multi::thrust::cuda::universal_array<complex, 3> fw_gpu(extents(in_gpu), 99.);

//      // multi::cuda::managed::array<complex, 3> const in_gpu(x, 10.);
//      // multi::cuda::managed::array<complex, 3> fw_gpu(extents(in_gpu), 99.);
//      {
//      //  boost::timer::auto_cpu_timer t; //  0.208237s wall, 0.200000s user + 0.010000s system = 0.210000s CPU (100.8%)
//          multi::cufft::dft({true, true}, in_gpu, fw_gpu, multi::cufft::forward);
//      //  BOOST_TEST( fw_gpu[8][9][10].operator complex() != 99. );
//      }
//      {
//      //  boost::timer::auto_cpu_timer t; //  0.208237s wall, 0.200000s user + 0.010000s system = 0.210000s CPU (100.8%)
//          multi::cufft::dft({true, true}, in_gpu, fw_gpu, multi::cufft::forward);
//      //  BOOST_TEST( fw_gpu[8][9][10].operator complex() != 99. );
//      }
//  }
// #endif
// }

// #if 0

// BOOST_AUTO_TEST_CASE(cufft_combinations, *utf::tolerance(0.00001)){

//  auto const in = []{
//      multi::array<complex, 4> ret({32, 90, 98, 96});
//      std::generate(ret.data_elements(), ret.data_elements() + ret.num_elements(),
//          [](){return complex{std::rand()*1./RAND_MAX, std::rand()*1./RAND_MAX};}
//      );
//      return ret;
//  }();
//  std::clog<<"memory size "<< in.num_elements()*sizeof(complex)/1e6 <<" MB\n";

//  multi::thrust::cuda::universal_array<complex, 4> const in_gpu = in;
//  multi::thrust::cuda::universal_array<complex, 4> const in_mng = in;

//  using std::clog;
//  for(auto c : std::vector<std::array<bool, 4>>{
//      {false, true , true , true },
//      {false, true , true , false},
//      {true , false, false, false},
//      {true , true , false, false},
//      {false, false, true , false},
//      {false, false, false, false},
//  }){
//      std::clog<<"case "; copy(begin(c), end(c), std::ostream_iterator<bool>{std::clog,", "}); std::clog<<std::endl;
//      multi::array<complex, 4> out = in;
//      multi::array<complex, 4> in_rw = in;
//      [&, _ = watch{"cpu_opl "}]{
//          multi::fftw::dft_forward(c, in, out);
//      }();
//      [&, _ = watch{"cpu_ipl "}]{
//          multi::fftw::dft(c, in_rw, multi::fftw::forward);
//      //  BOOST_TEST( abs( static_cast<multi::complex<double>>(in_rw[5][4][3][1]) - multi::complex<double>(out[5][4][3][1]) ) == 0. );
//      }();
//      {
//          multi::array<complex, 4> in_rw2 = in;
//          [&, _ = watch{"cpu_mov "}]{
//              multi::array<complex, 4> const out_mov = multi::fftw::dft_forward(c, std::move(in_rw2));
//          //  what(out_mov);
//          //  BOOST_TEST( abs( static_cast<multi::complex<double>>(out_mov[5][4][3][1]) - multi::complex<double>(out[5][4][3][1]) ) == 0. );
//              BOOST_REQUIRE( is_empty(in_rw2) );
//              BOOST_REQUIRE( extents(out_mov) == extents(in) );
//          }();
//      }

//      [&, _ = watch{"cpu_new "}]{
//          auto const out_cpy = multi::fftw::dft_forward(c, in);
//          BOOST_TEST( abs( static_cast<std::complex<double>>(out_cpy[5][4][3][1]) - std::complex<double>(out[5][4][3][1]) ) == 0. );
//      }();
//      multi::thrust::cuda::array<complex, 4> out_gpu(extents(in_gpu));
//      [&, _ = watch{"gpu_opl "}]{
//          multi::cufft::dft(c, in_gpu   , out_gpu, multi::cufft::forward);
//          BOOST_TEST( abs( static_cast<complex>(out_gpu[5][4][3][1]) - out[5][4][3][1] ) == 0. );
//      }();
//      {
//          multi::thrust::cuda::array<complex, 4> in_rw_gpu = in_gpu;
//          [&, _ = watch{"gpu_ipl "}]{
//              multi::cufft::dft(c, in_rw_gpu, multi::cufft::forward);
//              BOOST_TEST( abs( static_cast<complex>(in_rw_gpu[5][4][3][1]) - out[5][4][3][1] ) == 0. );
//          }();
//      }
//      {
//          multi::thrust::cuda::array<complex, 4> in_rw_gpu = in_gpu;
//          [&, _ = watch{"gpu_mov "}]{
//              multi::thrust::cuda::array<complex, 4> const out_mov = multi::cufft::dft_forward(c, std::move(in_rw_gpu));
//          //  BOOST_REQUIRE( in_rw_gpu.empty() );
//          //  BOOST_TEST( abs( static_cast<complex>(out_mov[5][4][3][1]) - out[5][4][3][1] ) == 0. );
//          }();
//      }
//      {
//          multi::thrust::cuda::array<complex, 4> in_rw_gpu = in_gpu;
//          [&, _ = watch{"gpu_mov "}]{
//              multi::thrust::cuda::array<complex, 4> out_mov = std::move(in_rw_gpu);
//              multi::cufft::dft(c, out_mov, multi::cufft::forward);
//          //  BOOST_REQUIRE( in_rw_gpu.empty() );
//          //  BOOST_TEST( abs( static_cast<complex>(out_mov[5][4][3][1]) - out[5][4][3][1] ) == 0. );
//          }();
//      }
//      cudaDeviceSynchronize();
//      [&, _ = watch{"gpu_new "}]{
//          multi::thrust::cuda::array<complex, 4> const out_cpy = multi::cufft::dft(c, in_gpu, multi::cufft::forward);
//      }();
//      multi::thrust::cuda::universal_array<complex, 4> out_mng(extents(in_mng));
//      [&, _ = watch{"mng_cld "}]{
//          multi::cufft::dft(c, in_mng, out_mng, multi::cufft::forward);
//          BOOST_TEST( abs( out_mng[5][4][3][1] - out[5][4][3][1] ) == 0. );
//      }();
//      [&, _ = watch{"mng_hot "}]{
//          multi::cufft::dft(c, in_mng   , out_mng, multi::cufft::forward);
//          BOOST_TEST( abs( out_mng[5][4][3][1] - out[5][4][3][1] ) == 0. );
//      }();
//      [&, _ = watch{"mng_new "}]{
//          auto const out_mng = multi::cufft::dft(c, in_mng, multi::cufft::forward);
//          BOOST_TEST( abs( out_mng[5][4][3][1] - out[5][4][3][1] ) == 0. );
//      }();
//  }
//  // std::clog<<"cache size "
//  //  << multi::cufft::plan::cache<1>().size() <<' '
//  //  << multi::cufft::plan::cache<2>().size() <<' '
//  //  << multi::cufft::plan::cache<3>().size() <<' '
//  //  << multi::cufft::plan::cache<4>().size() <<' '
//  // <<std::endl;
// }

// BOOST_AUTO_TEST_CASE(cufft_many_3D, *utf::tolerance(0.00001) ){

//  auto const in_cpu = []{
//      multi::array<complex, 4> ret({45, 18, 32, 16});
//      std::generate(
//          ret.data_elements(), ret.data_elements() + ret.num_elements(),
//          [](){return complex{std::rand()*1./RAND_MAX, std::rand()*1./RAND_MAX};}
//      );
//      return ret;
//  }();

//  multi::thrust::cuda::array<complex, 4> const in = in_cpu;
//  multi::thrust::cuda::array<complex, 4>       out(extents(in));

// #if 0
//  multi::cufft::many_dft(begin(unrotated(in)), end(unrotated(in)), begin(unrotated(out)), +1);

//  multi::array<complex, 4> out_cpu(extents(in));
//  multi::fft::many_dft(begin(unrotated(in_cpu)), end(unrotated(in_cpu)), begin(unrotated(out_cpu)), +1);

//  BOOST_TEST( imag( static_cast<complex>(out[5][4][3][2]) - out_cpu[5][4][3][2]) == 0. );
// #endif
// }

// #if 0
// BOOST_AUTO_TEST_CASE(cufft_4D, *utf::tolerance(0.00001) ){
//  auto const in = []{
//      multi::array<complex, 3> ret({10, 10, 10});
//      std::generate(ret.data_elements(), ret.data_elements() + ret.num_elements(),
//          [](){return complex{std::rand()*1./RAND_MAX, std::rand()*1./RAND_MAX};}
//      );
//      return ret;
//  }();

//  multi::array<complex, 3> out(extents(in));
// //  multi::fftw::dft({true, false, true}, in, out, multi::fftw::forward);
//  multi::fftw::many_dft(begin(in.rotated()), end(in.rotated()), begin(out.rotated()), multi::fftw::forward);

//  multi::thrust::cuda::array<complex, 3> in_gpu = in;
//  multi::thrust::cuda::array<complex, 3> out_gpu(extents(in));

// //  multi::cufft::dft({true, false, true}, in_gpu, out_gpu, multi::fft::forward);//multi::cufft::forward);
//  // multi::cufft::many_dft(begin(in_gpu.rotated()), end(in_gpu.rotated()), begin( out_gpu.rotated() ), multi::fftw::forward);
//  // BOOST_TEST( ( static_cast<complex>(out_gpu[5][4][3]) - out[5][4][3]).imag() == 0. );
// }
// #endif
// #endif

// #endif
