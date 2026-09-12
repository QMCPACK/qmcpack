// Copyright 2026 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/multi/adaptors/fftw.hpp>
#include <boost/multi/adaptors/thrust.hpp>
#include <boost/multi/adaptors/vkfft.hpp>
#include <boost/multi/array.hpp>

#include <boost/core/lightweight_test.hpp>

#include <thrust/complex.h>

#include <iostream>
#include <random>

namespace multi = boost::multi;
using complex   = thrust::complex<double>;

template<>
constexpr bool multi::force_element_trivial_default_construction<thrust::complex<double>> = true;

namespace {
auto random_array(multi::extents_t<2> exts) {
	multi::array<complex, 2>               ret(exts);
	std::mt19937                           gen(42);  // NOLINT(cert-msc32-c,cert-msc51-cpp) reproducible
	std::uniform_real_distribution<double> dist(-1.0, 1.0);
	std::generate(
		reinterpret_cast<double*>(ret.data_elements()),
		reinterpret_cast<double*>(ret.data_elements() + ret.num_elements()),
		[&] { return dist(gen); }
	);
	return ret;
}
}  // namespace

auto main() -> int {  // NOLINT(bugprone-exception-escape)
	complex const I{0.0, 1.0};  // NOLINT(readability-identifier-length)

	// 2D full C2C forward, out-of-place, cross-check against FFTW
	{
		auto const in_cpu = multi::array<complex, 2>{
			{ 1.0 + 2.0 * I,  9.0 - 1.0 * I,  2.0 + 4.0 * I},
			{ 3.0 + 3.0 * I,  7.0 - 4.0 * I,  1.0 + 9.0 * I},
			{ 4.0 + 1.0 * I,  5.0 + 3.0 * I,  2.0 + 4.0 * I},
			{ 3.0 - 1.0 * I,  8.0 + 7.0 * I,  2.0 + 1.0 * I},
			{31.0 - 1.0 * I, 18.0 + 7.0 * I, 2.0 + 10.0 * I},
		};

		auto fw_cpu = multi::array<complex, 2>(in_cpu.extents());
		multi::fftw::dft_forward({true, true}, in_cpu, fw_cpu);

		auto const in_gpu = multi::thrust::cuda::array<complex, 2>{in_cpu};
		auto       fw_gpu = multi::thrust::cuda::array<complex, 2>(in_gpu.extents());

		multi::vkfft::plan<2>({true, true}, in_gpu.layout(), fw_gpu.layout())
			.execute(in_gpu.base(), fw_gpu.base(), multi::vkfft::forward);

		BOOST_TEST(thrust::abs(complex(fw_gpu[3][2]) - fw_cpu[3][2]) < 1.0e-8);
		BOOST_TEST(thrust::abs(complex(fw_gpu[1][1]) - fw_cpu[1][1]) < 1.0e-8);
	}

	// 2D partial-axis masks (omitDimension)
	{
		auto const in_cpu = random_array({10, 20});

		for(auto which : std::vector<std::array<bool, 2>>{
				{ true,  true},
				{ true, false},
				{false,  true},
        }) {
			auto fw_cpu = multi::array<complex, 2>(in_cpu.extents());
			multi::fftw::dft(which, in_cpu, fw_cpu, multi::fftw::forward);

			auto const in_gpu = multi::thrust::cuda::array<complex, 2>{in_cpu};
			auto       fw_gpu = multi::thrust::cuda::array<complex, 2>(in_gpu.extents());

			multi::vkfft::dft(which, in_gpu, fw_gpu, multi::vkfft::forward);

			BOOST_TEST(thrust::abs(complex(fw_gpu[2][1]) - fw_cpu[2][1]) < 1.0e-8);
		}
	}

	// 3D full C2C forward
	{
		multi::array<complex, 3>               in_cpu({6, 8, 10});
		std::mt19937                           gen(7);  // NOLINT(cert-msc32-c,cert-msc51-cpp)
		std::uniform_real_distribution<double> dist(-1.0, 1.0);
		std::generate(in_cpu.elements().begin(), in_cpu.elements().end(), [&] { return complex{dist(gen), dist(gen)}; });

		auto fw_cpu = multi::array<complex, 3>(in_cpu.extents());
		multi::fftw::dft({true, true, true}, in_cpu, fw_cpu, multi::fftw::forward);

		auto const in_gpu = multi::thrust::cuda::array<complex, 3>{in_cpu};
		auto       fw_gpu = multi::thrust::cuda::array<complex, 3>(in_gpu.extents());

		multi::vkfft::dft({true, true, true}, in_gpu, fw_gpu, multi::vkfft::forward);

		BOOST_TEST(thrust::abs(complex(fw_gpu[3][2][1]) - fw_cpu[3][2][1]) < 1.0e-8);
	}

	// 3D = one batched axis + two transformed axes, for every placement of the
	// batch axis (VkFFT `omitDimension`).  Reference: FFTW with the same mask.
	{
		multi::array<complex, 3>               in_cpu({5, 12, 9});
		std::mt19937                           gen(11);  // NOLINT(cert-msc32-c,cert-msc51-cpp)
		std::uniform_real_distribution<double> dist(-1.0, 1.0);
		std::generate(in_cpu.elements().begin(), in_cpu.elements().end(), [&] { return complex{dist(gen), dist(gen)}; });

		auto const in_gpu = multi::thrust::cuda::array<complex, 3>{in_cpu};

		for(auto which : std::vector<std::array<bool, 3>>{
				{false,  true,  true},  // batch = outer axis
				{ true, false,  true},  // batch = middle axis
				{ true,  true, false},  // batch = inner axis (VkFFT "innermost batching")
        }) {
			auto fw_cpu = multi::array<complex, 3>(in_cpu.extents());
			multi::fftw::dft(which, in_cpu, fw_cpu, multi::fftw::forward);

			auto fw_gpu = multi::thrust::cuda::array<complex, 3>(in_gpu.extents());
			multi::vkfft::dft(which, in_gpu, fw_gpu, multi::vkfft::forward);

			multi::array<complex, 3> const fw_host = fw_gpu;

			double max_err = 0.0;
			auto   it_cpu  = fw_cpu.elements().begin();
			for(auto const& g : fw_host.elements()) { max_err = std::max(max_err, static_cast<double>(thrust::abs(g - *it_cpu++))); }
			BOOST_TEST(max_err < 1.0e-8);
		}
	}

	// 4D: two batched axes (outer + inner) + two transformed axes in the middle:
	// {false, true, true, false}
	{
		multi::array<complex, 4>               in_cpu({3, 10, 14, 4});
		std::mt19937                           gen(22);  // NOLINT(cert-msc32-c,cert-msc51-cpp)
		std::uniform_real_distribution<double> dist(-1.0, 1.0);
		std::generate(in_cpu.elements().begin(), in_cpu.elements().end(), [&] { return complex{dist(gen), dist(gen)}; });

		auto fw_cpu = multi::array<complex, 4>(in_cpu.extents());
		multi::fftw::dft({false, true, true, false}, in_cpu, fw_cpu, multi::fftw::forward);

		auto const in_gpu = multi::thrust::cuda::array<complex, 4>{in_cpu};
		auto       fw_gpu = multi::thrust::cuda::array<complex, 4>(in_gpu.extents());

		multi::vkfft::dft({false, true, true, false}, in_gpu, fw_gpu, multi::vkfft::forward);

		BOOST_TEST(thrust::abs(complex(fw_gpu[2][3][5][1]) - fw_cpu[2][3][5][1]) < 1.0e-8);
		BOOST_TEST(thrust::abs(complex(fw_gpu[0][9][0][3]) - fw_cpu[0][9][0][3]) < 1.0e-8);
	}

	// non-C-order layout: the contiguous axis is NOT the last one.
	// `.rotated()` of a contiguous X*Y*Z array has strides (Z, 1, Y*Z) -- the unit
	// stride is the middle axis.  The adaptor must sort axes by stride, not assume
	// axis D-1 is the fast one.  in and out share this rotated layout.
	{
		multi::array<complex, 3>               base_cpu({4, 5, 6});
		std::mt19937                           gen(33);  // NOLINT(cert-msc32-c,cert-msc51-cpp)
		std::uniform_real_distribution<double> dist(-1.0, 1.0);
		std::generate(base_cpu.elements().begin(), base_cpu.elements().end(), [&] { return complex{dist(gen), dist(gen)}; });

		auto fw_cpu = multi::array<complex, 3>(base_cpu.rotated().extents());  // contiguous 5x6x4
		multi::fftw::dft({true, true, true}, base_cpu.rotated(), fw_cpu, multi::fftw::forward);

		auto const base_gpu = multi::thrust::cuda::array<complex, 3>{base_cpu};
		auto       out_gpu  = multi::thrust::cuda::array<complex, 3>(base_gpu.extents());  // contiguous 4x5x6

		multi::vkfft::dft({true, true, true}, base_gpu.rotated(), out_gpu.rotated(), multi::vkfft::forward);

		multi::array<complex, 3> const out_host = out_gpu;
		double                         max_err  = 0.0;
		auto                           it_cpu   = fw_cpu.elements().begin();
		for(auto const& g : out_host.rotated().elements()) { max_err = std::max(max_err, static_cast<double>(thrust::abs(g - *it_cpu++))); }
		BOOST_TEST(max_err < 1.0e-8);
	}

	// in-place 2D
	{
		auto const in_cpu = random_array({8, 16});

		auto fw_cpu = in_cpu;
		multi::fftw::dft({true, true}, fw_cpu, multi::fftw::forward);

		auto fw_gpu = multi::thrust::cuda::array<complex, 2>{in_cpu};
		multi::vkfft::plan<2>({true, true}, fw_gpu.layout(), fw_gpu.layout())
			.execute(fw_gpu.base(), fw_gpu.base(), multi::vkfft::forward);

		BOOST_TEST(thrust::abs(complex(fw_gpu[2][1]) - fw_cpu[2][1]) < 1.0e-8);
	}

	// forward then backward round-trips to N * input (VkFFT inverse is unnormalized)
	{
		auto const in_cpu = random_array({12, 12});

		auto       roundtrip = multi::thrust::cuda::array<complex, 2>{in_cpu};
		auto       scratch   = multi::thrust::cuda::array<complex, 2>(roundtrip.extents());

		multi::vkfft::plan<2> const p({true, true}, roundtrip.layout(), scratch.layout());
		p.execute(roundtrip.base(), scratch.base(), multi::vkfft::forward);
		p.execute(scratch.base(), roundtrip.base(), multi::vkfft::backward);

		auto const n = static_cast<double>(in_cpu.num_elements());
		BOOST_TEST(thrust::abs(complex(roundtrip[5][7]) - n * in_cpu[5][7]) < 1.0e-6 * n);
	}

	// Exploratory: which `which` combinations does VkFFT reject / get wrong?
	//
	// Sweep every non-trivial mask in 1-D..4-D and compare against FFTW.  Any mask
	// where VkFFTAppend throws, or the result disagrees with FFTW, is reported.
	// (As of VkFFT v1.3.4 + VKFFT_MAX_FFT_DIMENSIONS=4 the interesting one is the
	// full 4-D transform {true,true,true,true}: FFTdim==4, which VkFFT's docs still
	// describe as "1, 2 or 3".)
	{
		auto sweep = [](auto exts, auto which) {
			constexpr auto D = std::tuple_size_v<decltype(which)>;
			multi::array<complex, D>                in_cpu(exts);
			std::mt19937                           gen(123);  // NOLINT(cert-msc32-c,cert-msc51-cpp)
			std::uniform_real_distribution<double> dist(-1.0, 1.0);
			std::generate(in_cpu.elements().begin(), in_cpu.elements().end(), [&] { return complex{dist(gen), dist(gen)}; });

			auto fw_cpu = multi::array<complex, D>(in_cpu.extents());
			multi::fftw::dft(which, in_cpu, fw_cpu, multi::fftw::forward);

			auto const in_gpu = multi::thrust::cuda::array<complex, D>{in_cpu};
			auto       fw_gpu = multi::thrust::cuda::array<complex, D>(in_gpu.extents());

			std::string label;
			for(bool b : which) { label += (b ? '1' : '0'); }

			try {
				multi::vkfft::dft(which, in_gpu, fw_gpu, multi::vkfft::forward);
			} catch(std::exception const& e) {
				std::cerr << "which={" << label << "}: VkFFT threw: " << e.what() << '\n';
				BOOST_TEST(false);
				return;
			}

			double max_err = 0.0;
			auto   it_cpu  = fw_cpu.elements().begin();
			for(auto const& g : fw_gpu.elements()) {
				max_err = std::max(max_err, static_cast<double>(thrust::abs(complex(g) - *it_cpu++)));
			}
			std::cerr << "which={" << label << "}: max|vkfft-fftw| = " << max_err << '\n';
			BOOST_TEST(max_err < 1.0e-7);
		};

		sweep(multi::extents_t<4>({3, 4, 5, 6}), std::array<bool, 4>{true, true, true, true});
		sweep(multi::extents_t<3>({4, 5, 6}), std::array<bool, 3>{true, false, true});
		sweep(multi::extents_t<4>({3, 4, 5, 6}), std::array<bool, 4>{true, false, false, true});
	}

	// The ONE `which` VkFFT (and this adaptor) refuses: the all-false mask.
	// VkFFT itself returns VKFFT_ERROR_UNSUPPORTED_FFT_OMIT (firstAxis > lastAxis)
	// when every dimension is omitted; the adaptor catches it earlier and throws
	// "vkfft: no transformed axis is not supported".  (FFTW instead treats an
	// all-false mask as a plain copy.)
	{
		auto const in_cpu = random_array({4, 4});
		auto const in_gpu = multi::thrust::cuda::array<complex, 2>{in_cpu};
		auto       fw_gpu = multi::thrust::cuda::array<complex, 2>(in_gpu.extents());

		bool threw = false;
		try {
			multi::vkfft::dft({false, false}, in_gpu, fw_gpu, multi::vkfft::forward);
		} catch(std::exception const& e) {
			threw = true;
			std::cerr << "which={00} (expected): " << e.what() << '\n';
		}
		BOOST_TEST(threw);
	}

	return boost::report_errors();
}
