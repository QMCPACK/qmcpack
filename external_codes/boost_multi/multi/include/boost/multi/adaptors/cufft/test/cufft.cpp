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
inline
#if defined(_MSC_VER) && !defined(__clang__)
__forceinline
#else
__attribute__((always_inline))
#endif
void DoNotOptimize(T const& value) {  // NOLINT(readability-identifier-naming) consistency with Google benchmark
#if defined(_MSC_VER) && !defined(__clang__)
	_ReadWriteBarrier(); (void)value;
#else
	asm volatile("" : "+m"(const_cast<T&>(value)));  // NOLINT(hicpp-no-assembler,cppcoreguidelines-pro-type-const-cast) hack
#endif
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

auto main() -> int {
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
			norm_t{}, 0.0, std::plus<>{}
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
}

