// Copyright 2021-2024 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/core/lightweight_test.hpp>

#include <boost/multi/adaptors/thrust.hpp>
#include <boost/multi/adaptors/thrust/managed_allocator.hpp>
#include <boost/multi/array.hpp>

#include <thrust/complex.h>
#include <thrust/device_allocator.h>
#include <thrust/device_ptr.h>
#include <thrust/system/cuda/memory.h>
#include <thrust/uninitialized_copy.h>
#include <thrust/universal_allocator.h>

#include <boost/timer/timer.hpp>

#include <numeric>

namespace multi = boost::multi;

#ifdef __NVCC__
template<>
inline constexpr bool ::boost::multi::force_element_trivial_default_construction<::std::complex<double>> = true;
template<>
inline constexpr bool ::boost::multi::force_element_trivial_default_construction<::std::complex<float>> = true;
template<>
inline constexpr bool ::boost::multi::force_element_trivial_default_construction<::thrust::complex<double>> = true;
template<>
inline constexpr bool ::boost::multi::force_element_trivial_default_construction<::thrust::complex<float>> = true;
#else  // vvv nvcc (12.1?) doesn't support this kind of customization: "error: expected initializer before ‘<’"
template<class T>
inline constexpr bool ::boost::multi::force_element_trivial_default_construction<::std::complex<T>> = std::is_trivially_default_constructible<T>::value;
template<class T>
inline constexpr bool ::boost::multi::force_element_trivial_default_construction<::thrust::complex<T>> = std::is_trivially_default_constructible<T>::value;
#endif

namespace {

template<class T> using test_allocator =
	//  multi::thrust::cuda::managed_allocator<T>
	thrust::cuda::allocator<T>;
}

// using types_list = boost::mpl::list<
//  // char,
//  double,
//  // std::complex<double>,
//  thrust::complex<double>
// >;

#define BOOST_AUTO_TEST_CASE(CasenamE) /**/

auto main() -> int {  // NOLINT(readability-function-cognitive-complexity,bugprone-exception-escape)

	BOOST_AUTO_TEST_CASE(cuda_universal_empty) {
		using complex = thrust::complex<double>;
		multi::array<complex, 2, thrust::cuda::universal_allocator<complex>> A;
		multi::array<complex, 2, thrust::cuda::universal_allocator<complex>> B = A;
		BOOST_TEST( A.is_empty() );
		BOOST_TEST( B.is_empty() );
		BOOST_TEST( A == B );
	}

	BOOST_AUTO_TEST_CASE(cuda_allocators) {

		multi::array<double, 1, thrust::cuda::allocator<double>> A1(200, 0.0);

		BOOST_TEST( size(A1) == 200 );
		A1[100] = 1.0;

		multi::array<double, 1, thrust::cuda::allocator<double>> const B1(200, 2.0);
		BOOST_TEST( B1[10] == 2.0 );

		A1[10] = B1[10];
		BOOST_TEST( A1[10] == 2.0 );
	}

	BOOST_AUTO_TEST_CASE(cuda_1d_initlist) {

		multi::array<double, 1, thrust::device_allocator<double>> A1 = {1.0, 2.0, 3.0};
		BOOST_TEST( A1.size() == 3 );

		// BOOST_TEST( size(A1) == 200 );
		// A1[100] = 1.0;

		// multi::array<double, 1, thrust::cuda::allocator<double>> const B1(200, 2.0);
		// BOOST_TEST( B1[10] == 2.0 );

		// A1[10] = B1[10];
		// BOOST_TEST( A1[10] == 2.0 );

		{
			multi::array<int, 1, thrust::device_allocator<int>> A = {1, 2, 3};
			multi::array<int, 1, thrust::device_allocator<int>> B(3, 0);

			BOOST_TEST( A[1] == 2 );

			thrust::transform(
				A.begin(), A.end(),
				B.begin(),
				[] __host__ __device__(int elem) {
					return elem * 2;
				}  // *2.0;}
			);

			BOOST_TEST( B[1] == 4 );
		}

		{
			multi::array<double, 1, thrust::device_allocator<double>> A = {1.0, 2.0, 3.0};
			multi::array<double, 1, thrust::device_allocator<double>> B(3);

			// // for(int i = 0; i != A.size(); ++i) { B[i] = A[i]*2.0; }
			// // for(auto i : A.extension()) { B[i] = A[i]*2.0; }

			thrust::transform(
				A.begin(), A.end(),
				B.begin(),
				[] __host__ __device__(double const& elem) { return elem * 2.0; }
			);

			BOOST_TEST( B[1] == 4.0 );
		}
	}

	BOOST_AUTO_TEST_CASE(test_univ_alloc) {
		multi::array<double, 2, thrust::cuda::universal_allocator<double>> Dev({128, 128});
		*raw_pointer_cast(Dev.base()) = 99.0;
	}

	BOOST_AUTO_TEST_CASE(mtc_universal_array) {
		multi::thrust::cuda::universal_array<double, 2> Dev({128, 128});
		*raw_pointer_cast(Dev.base()) = 99.0;
	}

	BOOST_AUTO_TEST_CASE(mtc_universal_coloncolon_array) {
		multi::thrust::cuda::universal::array<double, 2> Dev({128, 128});
		*raw_pointer_cast(Dev.base()) = 99.0;
	}

	BOOST_AUTO_TEST_CASE(test_alloc) {
		multi::array<double, 2, thrust::cuda::allocator<double>> Dev({128, 128});
		// *raw_pointer_cast(Dev.base()) = 99.0;  // segmentation fault (correct behavior)
	}

#ifdef NDEBUG
	BOOST_AUTO_TEST_CASE(thrust_copy_1D_issue123_double) {  // BOOST_AUTO_TEST_CASE(fdfdfdsfds) { using T = char;
		using T = double;

		static_assert(multi::is_trivially_default_constructible<T>{});
		static_assert(std::is_trivially_copy_constructible<T>{});
		static_assert(std::is_trivially_assignable<T&, T>{});

		multi::array<T, 1, test_allocator<T>> Devc(multi::extensions_t<1>{10240 * 10240});
		multi::array<T, 1, test_allocator<T>> Dev2(multi::extensions_t<1>{10240 * 10240});
		multi::array<T, 1>                    Host(multi::extensions_t<1>{10240 * 10240});
		std::iota(Host.elements().begin(), Host.elements().end(), 12.);
		multi::array<T, 1> Hos2(multi::extensions_t<1>{10240 * 10240});

		std::cout << "| 1D `" << typeid(T).name() << "` total data size: " << Host.num_elements() * sizeof(T) / 1073741824. << " GB | speed |\n|---|---|" << std::endl;
		{
			Devc = Host;
		}
		{
			boost::timer::auto_cpu_timer t{""};
			Devc = Host;
			cudaDeviceSynchronize();
			std::cout << "| contiguous host -> devc | " << Host.num_elements() * sizeof(T) / (t.elapsed().wall / 1e9) / 1073741824. << "GB/sec |" << std::endl;
		}
		{
			boost::timer::auto_cpu_timer t{""};
			Devc.sliced(0, 10240 * 10240 / 2) = Host.sliced(0, 10240 * 10240 / 2);
			cudaDeviceSynchronize();
			std::cout << "| sliced     host -> devc | " << Host.sliced(0, 10240 * 10240 / 2).num_elements() * sizeof(T) / (t.elapsed().wall / 1e9) / 1073741824. << "GB/sec |" << std::endl;
		}
		{
			boost::timer::auto_cpu_timer t{""};
			Devc.strided(2) = Host.strided(2);
			cudaDeviceSynchronize();
			std::cout << "| strided    host -> devc | " << Host.strided(2).num_elements() * sizeof(T) / (t.elapsed().wall / 1e9) / 1073741824. << "GB/sec |" << std::endl;
		}
		{
			boost::timer::auto_cpu_timer t{""};
			Hos2 = Devc;
			cudaDeviceSynchronize();
			std::cout << "| contiguous devc -> host | " << Host.num_elements() * sizeof(T) / (t.elapsed().wall / 1e9) / 1073741824. << "GB/sec |" << std::endl;
			// BOOST_TEST( Hos2 == Host );
		}
		{
			boost::timer::auto_cpu_timer t{""};
			Hos2.sliced(0, 10240 * 10240 / 2) = Devc.sliced(0, 10240 * 10240 / 2);
			cudaDeviceSynchronize();
			std::cout << "| sliced     devc -> host | " << Host.sliced(0, 10240 * 10240 / 2).num_elements() * sizeof(T) / (t.elapsed().wall / 1e9) / 1073741824. << "GB/sec |" << std::endl;
			// BOOST_TEST( Hos2 == Host );
		}
		{
			boost::timer::auto_cpu_timer t{""};
			Hos2.strided(2) = Devc.strided(2);
			cudaDeviceSynchronize();
			std::cout << "| strided    devc -> host | " << Host.strided(2).num_elements() * sizeof(T) / (t.elapsed().wall / 1e9) / 1073741824. << "GB/sec |" << std::endl;
			// BOOST_TEST( Hos2 == Host );
		}
		{
			boost::timer::auto_cpu_timer t{""};
			Dev2 = Devc;
			cudaDeviceSynchronize();
			std::cout << "| contiguous devc -> devc | " << Dev2.num_elements() * sizeof(T) / (t.elapsed().wall / 1e9) / 1073741824. << "GB/sec |" << std::endl;
			// BOOST_TEST( Dev2 == Devc );
		}
		{
			boost::timer::auto_cpu_timer t{""};
			auto                         Dev3 = Devc;
			cudaDeviceSynchronize();
			std::cout << "| copy_ctr   devc -> devc | " << Devc.num_elements() * sizeof(T) / (t.elapsed().wall / 1e9) / 1073741824. << "GB/sec |" << std::endl;
			// BOOST_TEST( Dev3 == Devc );
		}
		{
			boost::timer::auto_cpu_timer t{""};
			cudaMemcpy(raw_pointer_cast(Dev2.data_elements()), raw_pointer_cast(Devc.data_elements()), Devc.num_elements() * sizeof(T), cudaMemcpyDeviceToDevice);
			cudaDeviceSynchronize();
			std::cout << "| cudaMemcpy devc -> devc | " << Dev2.num_elements() * sizeof(T) / (t.elapsed().wall / 1e9) / 1073741824. << "GB/sec |" << std::endl;
			// BOOST_TEST( Dev2 == Devc );
		}
		{
			boost::timer::auto_cpu_timer t{""};
			cudaMemcpy(raw_pointer_cast(Dev2.data_elements()), raw_pointer_cast(Devc.data_elements()), Devc.num_elements() * sizeof(T), cudaMemcpyDeviceToDevice);
			cudaDeviceSynchronize();
			std::cout << "| cudaMemcpy devc -> devc | " << Dev2.num_elements() * sizeof(T) / (t.elapsed().wall / 1e9) / 1073741824. << "GB/sec |" << std::endl;
			// BOOST_TEST( Dev2 == Devc );
		}
		{
			boost::timer::auto_cpu_timer t{""};
			Dev2.sliced(0, 10240 * 10240 / 2) = Devc.sliced(0, 10240 * 10240 / 2);  //  0.005292s
			cudaDeviceSynchronize();
			std::cout << "| sliced     devc -> devc | " << Dev2.sliced(0, 10240 * 10240 / 2).num_elements() * sizeof(T) / (t.elapsed().wall / 1e9) / 1073741824. << "GB/sec |" << std::endl;
			// BOOST_TEST( Dev2 == Devc );
		}
		{
			boost::timer::auto_cpu_timer t{""};
			Dev2.strided(2) = Devc.strided(2);  //  0.005292s
			cudaDeviceSynchronize();
			std::cout << "| strided    devc -> devc | " << Dev2.strided(2).num_elements() * sizeof(T) / (t.elapsed().wall / 1e9) / 1073741824. << "GB/sec |" << std::endl;
			// BOOST_TEST( Dev2 == Devc );
		}
		{
			boost::timer::auto_cpu_timer t{""};
			Hos2 = Host;
			cudaDeviceSynchronize();
			std::cout << "| contiguous host -> host | " << Hos2.num_elements() * sizeof(T) / (t.elapsed().wall / 1e9) / 1073741824. << "GB/sec |" << std::endl;
		}
		{
			boost::timer::auto_cpu_timer t{""};
			Hos2.sliced(0, 10240 * 10240 / 2) = Host.sliced(0, 10240 * 10240 / 2);  //  0.005292s
			cudaDeviceSynchronize();
			std::cout << "| sliced     host -> host | " << Hos2.sliced(0, 10240 * 10240 / 2).num_elements() * sizeof(T) / (t.elapsed().wall / 1e9) / 1073741824. << "GB/sec |" << std::endl;
		}
		{
			boost::timer::auto_cpu_timer t{""};
			Hos2.strided(2) = Host.strided(2);  //  0.005292s
			cudaDeviceSynchronize();
			std::cout << "| strided    host -> host | " << Hos2.strided(2).num_elements() * sizeof(T) / (t.elapsed().wall / 1e9) / 1073741824. << "GB/sec |" << std::endl;
		}
		std::cout << "   " << std::endl;
	}

	BOOST_AUTO_TEST_CASE(thrust_copy_1D_issue123_complex) {  // BOOST_AUTO_TEST_CASE(fdfdfdsfds) { using T = char;
		using T = thrust::complex<double>;

		static_assert(multi::is_trivially_default_constructible<T>{}, "!");
		static_assert(std::is_trivially_copy_constructible<T>{}, "!");
		static_assert(std::is_trivially_assignable<T&, T>{}, "!");

		multi::array<T, 1, test_allocator<T>> Devc(multi::extensions_t<1>{10240 * 10240});
		multi::array<T, 1, test_allocator<T>> Dev2(multi::extensions_t<1>{10240 * 10240});
		multi::array<T, 1>                    Host(multi::extensions_t<1>{10240 * 10240});
		std::iota(Host.elements().begin(), Host.elements().end(), 12.);
		multi::array<T, 1> Hos2(multi::extensions_t<1>{10240 * 10240});

		std::cout << "| 1D `" << typeid(T).name() << "` total data size: " << Host.num_elements() * sizeof(T) / 1073741824. << " GB | speed |\n|---|---|" << std::endl;
		{
			Devc = Host;
		}
		{
			boost::timer::auto_cpu_timer t{""};
			Devc = Host;
			cudaDeviceSynchronize();
			std::cout << "| contiguous host -> devc | " << Host.num_elements() * sizeof(T) / (t.elapsed().wall / 1e9) / 1073741824. << "GB/sec |" << std::endl;
		}
		{
			boost::timer::auto_cpu_timer t{""};
			Devc.sliced(0, 10240 * 10240 / 2) = Host.sliced(0, 10240 * 10240 / 2);
			cudaDeviceSynchronize();
			std::cout << "| sliced     host -> devc | " << Host.sliced(0, 10240 * 10240 / 2).num_elements() * sizeof(T) / (t.elapsed().wall / 1e9) / 1073741824. << "GB/sec |" << std::endl;
		}
		{
			boost::timer::auto_cpu_timer t{""};
			Devc.strided(2) = Host.strided(2);
			cudaDeviceSynchronize();
			std::cout << "| strided    host -> devc | " << Host.strided(2).num_elements() * sizeof(T) / (t.elapsed().wall / 1e9) / 1073741824. << "GB/sec |" << std::endl;
		}
		{
			boost::timer::auto_cpu_timer t{""};
			Hos2 = Devc;
			cudaDeviceSynchronize();
			std::cout << "| contiguous devc -> host | " << Host.num_elements() * sizeof(T) / (t.elapsed().wall / 1e9) / 1073741824. << "GB/sec |" << std::endl;
			// BOOST_TEST( Hos2 == Host );
		}
		{
			boost::timer::auto_cpu_timer t{""};
			Hos2.sliced(0, 10240 * 10240 / 2) = Devc.sliced(0, 10240 * 10240 / 2);
			cudaDeviceSynchronize();
			std::cout << "| sliced     devc -> host | " << Host.sliced(0, 10240 * 10240 / 2).num_elements() * sizeof(T) / (t.elapsed().wall / 1e9) / 1073741824. << "GB/sec |" << std::endl;
			// BOOST_TEST( Hos2 == Host );
		}
		{
			boost::timer::auto_cpu_timer t{""};
			Hos2.strided(2) = Devc.strided(2);
			cudaDeviceSynchronize();
			std::cout << "| strided    devc -> host | " << Host.strided(2).num_elements() * sizeof(T) / (t.elapsed().wall / 1e9) / 1073741824. << "GB/sec |" << std::endl;
			// BOOST_TEST( Hos2 == Host );
		}
		{
			boost::timer::auto_cpu_timer t{""};
			Dev2 = Devc;
			cudaDeviceSynchronize();
			std::cout << "| contiguous devc -> devc | " << Dev2.num_elements() * sizeof(T) / (t.elapsed().wall / 1e9) / 1073741824. << "GB/sec |" << std::endl;
			// BOOST_TEST( Dev2 == Devc );
		}
		{
			boost::timer::auto_cpu_timer t{""};
			auto                         Dev3 = Devc;
			cudaDeviceSynchronize();
			std::cout << "| copy_ctr   devc -> devc | " << Devc.num_elements() * sizeof(T) / (t.elapsed().wall / 1e9) / 1073741824. << "GB/sec |" << std::endl;
			// BOOST_TEST( Dev3 == Devc );
		}
		{
			boost::timer::auto_cpu_timer t{""};
			cudaMemcpy(raw_pointer_cast(Dev2.data_elements()), raw_pointer_cast(Devc.data_elements()), Devc.num_elements() * sizeof(T), cudaMemcpyDeviceToDevice);
			cudaDeviceSynchronize();
			std::cout << "| cudaMemcpy devc -> devc | " << Dev2.num_elements() * sizeof(T) / (t.elapsed().wall / 1e9) / 1073741824. << "GB/sec |" << std::endl;
			// BOOST_TEST( Dev2 == Devc );
		}
		{
			boost::timer::auto_cpu_timer t{""};
			cudaMemcpy(raw_pointer_cast(Dev2.data_elements()), raw_pointer_cast(Devc.data_elements()), Devc.num_elements() * sizeof(T), cudaMemcpyDeviceToDevice);
			cudaDeviceSynchronize();
			std::cout << "| cudaMemcpy devc -> devc | " << Dev2.num_elements() * sizeof(T) / (t.elapsed().wall / 1e9) / 1073741824. << "GB/sec |" << std::endl;
			// BOOST_TEST( Dev2 == Devc );
		}
		{
			boost::timer::auto_cpu_timer t{""};
			Dev2.sliced(0, 10240 * 10240 / 2) = Devc.sliced(0, 10240 * 10240 / 2);  //  0.005292s
			cudaDeviceSynchronize();
			std::cout << "| sliced     devc -> devc | " << Dev2.sliced(0, 10240 * 10240 / 2).num_elements() * sizeof(T) / (t.elapsed().wall / 1e9) / 1073741824. << "GB/sec |" << std::endl;
			// BOOST_TEST( Dev2 == Devc );
		}
		{
			boost::timer::auto_cpu_timer t{""};
			Dev2.strided(2) = Devc.strided(2);  //  0.005292s
			cudaDeviceSynchronize();
			std::cout << "| strided    devc -> devc | " << Dev2.strided(2).num_elements() * sizeof(T) / (t.elapsed().wall / 1e9) / 1073741824. << "GB/sec |" << std::endl;
			// BOOST_TEST( Dev2 == Devc );
		}
		{
			boost::timer::auto_cpu_timer t{""};
			Hos2 = Host;
			cudaDeviceSynchronize();
			std::cout << "| contiguous host -> host | " << Hos2.num_elements() * sizeof(T) / (t.elapsed().wall / 1e9) / 1073741824. << "GB/sec |" << std::endl;
		}
		{
			boost::timer::auto_cpu_timer t{""};
			Hos2.sliced(0, 10240 * 10240 / 2) = Host.sliced(0, 10240 * 10240 / 2);  //  0.005292s
			cudaDeviceSynchronize();
			std::cout << "| sliced     host -> host | " << Hos2.sliced(0, 10240 * 10240 / 2).num_elements() * sizeof(T) / (t.elapsed().wall / 1e9) / 1073741824. << "GB/sec |" << std::endl;
		}
		{
			boost::timer::auto_cpu_timer t{""};
			Hos2.strided(2) = Host.strided(2);  //  0.005292s
			cudaDeviceSynchronize();
			std::cout << "| strided    host -> host | " << Hos2.strided(2).num_elements() * sizeof(T) / (t.elapsed().wall / 1e9) / 1073741824. << "GB/sec |" << std::endl;
		}
		std::cout << "   " << std::endl;
	}

	BOOST_AUTO_TEST_CASE(thrust_cpugpu_2D_issue123_double) {
		using T = double;

		auto const exts = multi::extensions_t<2>({10240, 10240});

		std::cout << "| 2D `" << typeid(T).name() << "` max data size " << exts.num_elements() * sizeof(T) / 1073741824. << " GB | speed |\n|---|---|" << std::endl;

		multi::array<T, 2, test_allocator<T>> Devc(exts);
		multi::array<T, 2, test_allocator<T>> Dev2(exts);

		multi::array<T, 2> Host(exts);
		std::iota(Host.elements().begin(), Host.elements().end(), 12.);
		multi::array<T, 2> Hos2(exts);

		{
			Devc({0, 5120}, {0, 5120}) = Host({0, 5120}, {0, 5120});  // 0.002859s
		}
		{
			boost::timer::auto_cpu_timer t{""};
			Devc = Host;
			std::cout << "| contiguous host to devc | " << Host.num_elements() * sizeof(T) / (t.elapsed().wall / 1e9) / 1073741824. << "GB/sec |" << std::endl;
		}
		{
			boost::timer::auto_cpu_timer t{""};
			Devc.sliced(0, 5120) = Host.sliced(0, 5120);  //  0.005292s
			std::cout << "| sliced     host to devc | " << Host.sliced(0, 5120).num_elements() * sizeof(T) / (t.elapsed().wall / 1e9) / 1073741824. << "GB/sec |" << std::endl;
		}
		{
			boost::timer::auto_cpu_timer t{""};
			Devc({0, 5120}, {0, 5120}) = Host({0, 5120}, {0, 5120});  // 0.002859s
			std::cout << "| strided    host to devc | " << Host({0, 5120}, {0, 5120}).num_elements() * sizeof(T) / (t.elapsed().wall / 1e9) / 1073741824. << "GB/sec |" << std::endl;
		}
		{
			boost::timer::auto_cpu_timer t{""};
			Hos2 = Devc;
			std::cout << "| contiguous devc to host | " << Host.num_elements() * sizeof(T) / (t.elapsed().wall / 1e9) / 1073741824. << "GB/sec |" << std::endl;
			// BOOST_TEST(Hos2 == Host);
		}
		{
			boost::timer::auto_cpu_timer t{""};
			Hos2.sliced(0, 5120) = Devc.sliced(0, 5120);  //  0.005292s
			std::cout << "| sliced     devc to host | " << Host.sliced(0, 5120).num_elements() * sizeof(T) / (t.elapsed().wall / 1e9) / 1073741824. << "GB/sec |" << std::endl;
			// BOOST_TEST(Hos2 == Host);
		}
		{
			boost::timer::auto_cpu_timer t{""};
			Hos2({0, 5120}, {0, 5120}) = Devc({0, 5120}, {0, 5120});  // 0.002859s
			std::cout << "| strided    devc to host | " << Host({0, 5120}, {0, 5120}).num_elements() * sizeof(T) / (t.elapsed().wall / 1e9) / 1073741824. << "GB/sec |" << std::endl;
			// BOOST_TEST(Hos2 == Host);
		}
		{
			boost::timer::auto_cpu_timer t{""};
			Dev2 = Devc;
			cudaDeviceSynchronize();
			std::cout << "| contiguous devc to devc | " << Host.num_elements() * sizeof(T) / (t.elapsed().wall / 1e9) / 1073741824. << "GB/sec |" << std::endl;
			// BOOST_TEST( Dev2 == Devc );
		}
		{
			boost::timer::auto_cpu_timer t{""};
			auto                         Dev3 = Devc;
			cudaDeviceSynchronize();
			std::cout << "| copy_ctr   devc -> devc | " << Devc.num_elements() * sizeof(T) / (t.elapsed().wall / 1e9) / 1073741824. << "GB/sec |" << std::endl;
			// BOOST_TEST( Dev3 == Devc );
		}
		{
			boost::timer::auto_cpu_timer t{""};
			cudaMemcpy(raw_pointer_cast(Dev2.data_elements()), raw_pointer_cast(Devc.data_elements()), Devc.num_elements() * sizeof(T), cudaMemcpyDeviceToDevice);
			cudaDeviceSynchronize();
			std::cout << "| cudaMemcpy devc to devc | " << Host.num_elements() * sizeof(T) / (t.elapsed().wall / 1e9) / 1073741824. << "GB/sec |" << std::endl;
			// BOOST_TEST( Dev2 == Devc );
		}
		{
			boost::timer::auto_cpu_timer t{""};
			cudaMemcpy(raw_pointer_cast(Dev2.data_elements()), raw_pointer_cast(Devc.data_elements()), Devc.num_elements() * sizeof(T), cudaMemcpyDeviceToDevice);
			cudaDeviceSynchronize();
			std::cout << "| cudaMemcpy devc to devc | " << Host.num_elements() * sizeof(T) / (t.elapsed().wall / 1e9) / 1073741824. << "GB/sec |" << std::endl;
			// BOOST_TEST( Dev2 == Devc );
		}
		{
			boost::timer::auto_cpu_timer t{""};
			auto                         Dev3 = Devc;
			cudaDeviceSynchronize();
			std::cout << "| copy_ctr   devc to devc | " << Host.num_elements() * sizeof(T) / (t.elapsed().wall / 1e9) / 1073741824. << "GB/sec |" << std::endl;
			// BOOST_TEST( Dev2 == Devc );
		}
		{
			boost::timer::auto_cpu_timer t{""};
			Dev2.sliced(0, 5120) = Devc.sliced(0, 5120);  //  0.005292s
			cudaDeviceSynchronize();
			std::cout << "| sliced     devc to devc | " << Host.sliced(0, 5120).num_elements() * sizeof(T) / (t.elapsed().wall / 1e9) / 1073741824. << "GB/sec |" << std::endl;
			// BOOST_TEST( Dev2 == Devc );
		}
		{
			boost::timer::auto_cpu_timer t{""};
			Dev2({0, 5120}, {0, 5120}) = Devc({0, 5120}, {0, 5120});  // 0.002859s
			cudaDeviceSynchronize();
			std::cout << "| strided    devc to devc | " << Host({0, 5120}, {0, 5120}).num_elements() * sizeof(T) / (t.elapsed().wall / 1e9) / 1073741824. << "GB/sec |" << std::endl;
			// BOOST_TEST( Dev2 == Devc );
		}
		{
			boost::timer::auto_cpu_timer t{""};
			Hos2 = Host;
			std::cout << "| contiguous host to host | " << Host.num_elements() * sizeof(T) / (t.elapsed().wall / 1e9) / 1073741824. << "GB/sec |" << std::endl;
		}
		{
			boost::timer::auto_cpu_timer t{""};
			Hos2.sliced(0, 5120) = Host.sliced(0, 5120);  //  0.005292s
			std::cout << "| sliced     host to host | " << Host.sliced(0, 5120).num_elements() * sizeof(T) / (t.elapsed().wall / 1e9) / 1073741824. << "GB/sec |" << std::endl;
		}
		{
			boost::timer::auto_cpu_timer t{""};
			Hos2({0, 5120}, {0, 5120}) = Host({0, 5120}, {0, 5120});  // 0.002859s
			std::cout << "| strided    host to host | " << Host({0, 5120}, {0, 5120}).num_elements() * sizeof(T) / (t.elapsed().wall / 1e9) / 1073741824. << "GB/sec |" << std::endl;
		}
		std::cout << "  " << std::endl;
	}

	BOOST_AUTO_TEST_CASE(thrust_cpugpu_2D_issue123_complex) {
		using T = thrust::complex<double>;

		auto const exts = multi::extensions_t<2>({10240, 10240});

		std::cout << "| 2D `" << typeid(T).name() << "` max data size " << exts.num_elements() * sizeof(T) / 1073741824. << " GB | speed |\n|---|---|" << std::endl;

		multi::array<T, 2, test_allocator<T>> Devc(exts);
		multi::array<T, 2, test_allocator<T>> Dev2(exts);

		multi::array<T, 2> Host(exts);
		std::iota(Host.elements().begin(), Host.elements().end(), 12.);
		multi::array<T, 2> Hos2(exts);

		{
			Devc({0, 5120}, {0, 5120}) = Host({0, 5120}, {0, 5120});  // 0.002859s
		}
		{
			boost::timer::auto_cpu_timer t{""};
			Devc = Host;
			std::cout << "| contiguous host to devc | " << Host.num_elements() * sizeof(T) / (t.elapsed().wall / 1e9) / 1073741824. << "GB/sec |" << std::endl;
		}
		{
			boost::timer::auto_cpu_timer t{""};
			Devc.sliced(0, 5120) = Host.sliced(0, 5120);  //  0.005292s
			std::cout << "| sliced     host to devc | " << Host.sliced(0, 5120).num_elements() * sizeof(T) / (t.elapsed().wall / 1e9) / 1073741824. << "GB/sec |" << std::endl;
		}
		{
			boost::timer::auto_cpu_timer t{""};
			Devc({0, 5120}, {0, 5120}) = Host({0, 5120}, {0, 5120});  // 0.002859s
			std::cout << "| strided    host to devc | " << Host({0, 5120}, {0, 5120}).num_elements() * sizeof(T) / (t.elapsed().wall / 1e9) / 1073741824. << "GB/sec |" << std::endl;
		}
		{
			boost::timer::auto_cpu_timer t{""};
			Hos2 = Devc;
			std::cout << "| contiguous devc to host | " << Host.num_elements() * sizeof(T) / (t.elapsed().wall / 1e9) / 1073741824. << "GB/sec |" << std::endl;
			// BOOST_TEST(Hos2 == Host);
		}
		{
			boost::timer::auto_cpu_timer t{""};
			Hos2.sliced(0, 5120) = Devc.sliced(0, 5120);  //  0.005292s
			std::cout << "| sliced     devc to host | " << Host.sliced(0, 5120).num_elements() * sizeof(T) / (t.elapsed().wall / 1e9) / 1073741824. << "GB/sec |" << std::endl;
			// BOOST_TEST(Hos2 == Host);
		}
		{
			boost::timer::auto_cpu_timer t{""};
			Hos2({0, 5120}, {0, 5120}) = Devc({0, 5120}, {0, 5120});  // 0.002859s
			std::cout << "| strided    devc to host | " << Host({0, 5120}, {0, 5120}).num_elements() * sizeof(T) / (t.elapsed().wall / 1e9) / 1073741824. << "GB/sec |" << std::endl;
			// BOOST_TEST(Hos2 == Host);
		}
		{
			boost::timer::auto_cpu_timer t{""};
			Dev2 = Devc;
			cudaDeviceSynchronize();
			std::cout << "| contiguous devc to devc | " << Host.num_elements() * sizeof(T) / (t.elapsed().wall / 1e9) / 1073741824. << "GB/sec |" << std::endl;
			// BOOST_TEST( Dev2 == Devc );
		}
		{
			boost::timer::auto_cpu_timer t{""};
			auto                         Dev3 = Devc;
			cudaDeviceSynchronize();
			std::cout << "| copy_ctr   devc -> devc | " << Devc.num_elements() * sizeof(T) / (t.elapsed().wall / 1e9) / 1073741824. << "GB/sec |" << std::endl;
			// BOOST_TEST( Dev3 == Devc );
		}
		{
			boost::timer::auto_cpu_timer t{""};
			cudaMemcpy(raw_pointer_cast(Dev2.data_elements()), raw_pointer_cast(Devc.data_elements()), Devc.num_elements() * sizeof(T), cudaMemcpyDeviceToDevice);
			cudaDeviceSynchronize();
			std::cout << "| cudaMemcpy devc to devc | " << Host.num_elements() * sizeof(T) / (t.elapsed().wall / 1e9) / 1073741824. << "GB/sec |" << std::endl;
			// BOOST_TEST( Dev2 == Devc );
		}
		{
			boost::timer::auto_cpu_timer t{""};
			cudaMemcpy(raw_pointer_cast(Dev2.data_elements()), raw_pointer_cast(Devc.data_elements()), Devc.num_elements() * sizeof(T), cudaMemcpyDeviceToDevice);
			cudaDeviceSynchronize();
			std::cout << "| cudaMemcpy devc to devc | " << Host.num_elements() * sizeof(T) / (t.elapsed().wall / 1e9) / 1073741824. << "GB/sec |" << std::endl;
			// BOOST_TEST( Dev2 == Devc );
		}
		{
			boost::timer::auto_cpu_timer t{""};
			auto                         Dev3 = Devc;
			cudaDeviceSynchronize();
			std::cout << "| copy_ctr   devc to devc | " << Host.num_elements() * sizeof(T) / (t.elapsed().wall / 1e9) / 1073741824. << "GB/sec |" << std::endl;
			// BOOST_TEST( Dev2 == Devc );
		}
		{
			boost::timer::auto_cpu_timer t{""};
			Dev2.sliced(0, 5120) = Devc.sliced(0, 5120);  //  0.005292s
			cudaDeviceSynchronize();
			std::cout << "| sliced     devc to devc | " << Host.sliced(0, 5120).num_elements() * sizeof(T) / (t.elapsed().wall / 1e9) / 1073741824. << "GB/sec |" << std::endl;
			// BOOST_TEST( Dev2 == Devc );
		}
		{
			boost::timer::auto_cpu_timer t{""};
			Dev2({0, 5120}, {0, 5120}) = Devc({0, 5120}, {0, 5120});  // 0.002859s
			cudaDeviceSynchronize();
			std::cout << "| strided    devc to devc | " << Host({0, 5120}, {0, 5120}).num_elements() * sizeof(T) / (t.elapsed().wall / 1e9) / 1073741824. << "GB/sec |" << std::endl;
			// BOOST_TEST( Dev2 == Devc );
		}
		{
			boost::timer::auto_cpu_timer t{""};
			Hos2 = Host;
			std::cout << "| contiguous host to host | " << Host.num_elements() * sizeof(T) / (t.elapsed().wall / 1e9) / 1073741824. << "GB/sec |" << std::endl;
		}
		{
			boost::timer::auto_cpu_timer t{""};
			Hos2.sliced(0, 5120) = Host.sliced(0, 5120);  //  0.005292s
			std::cout << "| sliced     host to host | " << Host.sliced(0, 5120).num_elements() * sizeof(T) / (t.elapsed().wall / 1e9) / 1073741824. << "GB/sec |" << std::endl;
		}
		{
			boost::timer::auto_cpu_timer t{""};
			Hos2({0, 5120}, {0, 5120}) = Host({0, 5120}, {0, 5120});  // 0.002859s
			std::cout << "| strided    host to host | " << Host({0, 5120}, {0, 5120}).num_elements() * sizeof(T) / (t.elapsed().wall / 1e9) / 1073741824. << "GB/sec |" << std::endl;
		}
		std::cout << "  " << std::endl;
	}

	BOOST_AUTO_TEST_CASE(thrust_cpugpu_issue123_3D_double) {
		using T         = double;
		auto const exts = multi::extensions_t<3>({1024, 1024, 100});

		std::cout << "| 3D `" << typeid(T).name() << "` max data size " << exts.num_elements() * sizeof(T) / 1073741824. << " GB | speed |\n|---|---|" << std::endl;

		multi::array<T, 3, test_allocator<T>> Devc(exts);
		multi::array<T, 3, test_allocator<T>> Dev2(exts);
		multi::array<T, 3>                    Host(exts);
		std::iota(Host.elements().begin(), Host.elements().end(), 12.);
		multi::array<T, 3> Hos2(exts);

		{
			Devc({0, 512}, {0, 512}, {0, 512}) = Host({0, 512}, {0, 512}, {0, 512});  // 0.002859s
		}
		{
			boost::timer::auto_cpu_timer t{""};
			Devc = Host;
			std::cout << "| contiguous host to devc | " << Host.num_elements() * sizeof(T) / (t.elapsed().wall / 1e9) / 1073741824. << " GB/sec |" << std::endl;
		}
		{
			boost::timer::auto_cpu_timer t{""};
			Devc.sliced(0, 512) = Host.sliced(0, 512);  //  0.005292s
			std::cout << "| sliced     host to devc | " << Host.sliced(0, 512).num_elements() * sizeof(T) / (t.elapsed().wall / 1e9) / 1073741824. << " GB/sec |" << std::endl;
		}
		{
			boost::timer::auto_cpu_timer t{""};
			Devc({0, 512}, {0, 512}, {0, 512}) = Host({0, 512}, {0, 512}, {0, 512});  // 0.002859s
			std::cout << "| strided    host to devc | " << Host({0, 512}, {0, 512}, {0, 512}).num_elements() * sizeof(T) / (t.elapsed().wall / 1e9) / 1073741824. << "GB/sec |" << std::endl;
		}
		{
			boost::timer::auto_cpu_timer t{""};
			Hos2 = Devc;
			std::cout << "| contiguous devc to host | " << Host.num_elements() * sizeof(T) / (t.elapsed().wall / 1e9) / 1073741824. << "GB/sec |" << std::endl;
			// BOOST_TEST(Hos2 == Host);
		}
		{
			boost::timer::auto_cpu_timer t{""};
			Hos2.sliced(0, 512) = Devc.sliced(0, 512);  //  0.005292s
			std::cout << "| sliced     devc to host | " << Host.sliced(0, 512).num_elements() * sizeof(T) / (t.elapsed().wall / 1e9) / 1073741824. << "GB/sec |" << std::endl;
			// BOOST_TEST(Hos2 == Host);
		}
		{
			boost::timer::auto_cpu_timer t{""};
			Hos2({0, 512}, {0, 512}, {0, 512}) = Devc({0, 512}, {0, 512}, {0, 512});  // 0.002859s
			std::cout << "| strided    devc to host | " << Host({0, 512}, {0, 512}, {0, 512}).num_elements() * sizeof(T) / (t.elapsed().wall / 1e9) / 1073741824. << "GB/sec |" << std::endl;
			// BOOST_TEST(Hos2 == Host);
		}
		{
			boost::timer::auto_cpu_timer t{""};
			Dev2 = Devc;
			cudaDeviceSynchronize();
			std::cout << "| contiguous devc to devc | " << Dev2.num_elements() * sizeof(T) / (t.elapsed().wall / 1e9) / 1073741824. << "GB/sec |" << std::endl;
			// BOOST_TEST(Dev2 == Devc);
		}
		{
			boost::timer::auto_cpu_timer t{""};
			auto                         Dev3 = Devc;
			cudaDeviceSynchronize();
			std::cout << "| copy_ctr   devc -> devc | " << Devc.num_elements() * sizeof(T) / (t.elapsed().wall / 1e9) / 1073741824. << "GB/sec |" << std::endl;
			// BOOST_TEST( Dev3 == Devc );
		}
		{
			boost::timer::auto_cpu_timer t{""};
			cudaMemcpy(raw_pointer_cast(Dev2.data_elements()), raw_pointer_cast(Devc.data_elements()), Devc.num_elements() * sizeof(T), cudaMemcpyDeviceToDevice);
			cudaDeviceSynchronize();
			std::cout << "| cudaMemcpy devc -> devc | " << Dev2.num_elements() * sizeof(T) / (t.elapsed().wall / 1e9) / 1073741824. << "GB/sec |" << std::endl;
			// BOOST_TEST( Dev2 == Devc );
		}
		{
			boost::timer::auto_cpu_timer t{""};
			cudaMemcpy(raw_pointer_cast(Dev2.data_elements()), raw_pointer_cast(Devc.data_elements()), Devc.num_elements() * sizeof(T), cudaMemcpyDeviceToDevice);
			cudaDeviceSynchronize();
			std::cout << "| cudaMemcpy devc -> devc | " << Dev2.num_elements() * sizeof(T) / (t.elapsed().wall / 1e9) / 1073741824. << "GB/sec |" << std::endl;
			// BOOST_TEST( Dev2 == Devc );
		}
		{
			boost::timer::auto_cpu_timer t{""};
			Dev2.sliced(0, 512) = Devc.sliced(0, 512);  //  0.005292s
			cudaDeviceSynchronize();
			std::cout << "| sliced     devc to devc | " << Dev2.sliced(0, 512).num_elements() * sizeof(T) / (t.elapsed().wall / 1e9) / 1073741824. << "GB/sec |" << std::endl;
			// BOOST_TEST(Dev2 == Devc);
		}
		{
			boost::timer::auto_cpu_timer t{""};
			Dev2({0, 512}, {0, 512}, {0, 512}) = Devc({0, 512}, {0, 512}, {0, 512});  // 0.002859s
			cudaDeviceSynchronize();
			std::cout << "| strided    devc to devc | " << Dev2({0, 512}, {0, 512}, {0, 512}).num_elements() * sizeof(T) / (t.elapsed().wall / 1e9) / 1073741824. << "GB/sec |" << std::endl;
			// BOOST_TEST(Dev2 == Devc);
		}
		{
			boost::timer::auto_cpu_timer t{""};
			Hos2 = Host;
			std::cout << "| contiguous host to host | " << Hos2.num_elements() * sizeof(T) / (t.elapsed().wall / 1e9) / 1073741824. << "GB/sec |" << std::endl;
		}
		{
			boost::timer::auto_cpu_timer t{""};
			Hos2.sliced(0, 512) = Host.sliced(0, 512);  //  0.005292s
			std::cout << "| sliced     host to host | " << Hos2.sliced(0, 512).num_elements() * sizeof(T) / (t.elapsed().wall / 1e9) / 1073741824. << "GB/sec |" << std::endl;
		}
		{
			boost::timer::auto_cpu_timer t{""};
			Hos2({0, 512}, {0, 512}, {0, 512}) = Host({0, 512}, {0, 512}, {0, 512});  // 0.002859s
			std::cout << "| strided    host to host | " << Hos2({0, 512}, {0, 512}, {0, 512}).num_elements() * sizeof(T) / (t.elapsed().wall / 1e9) / 1073741824. << "GB/sec |" << std::endl;
		}
		std::cout << "   " << std::endl;
	}

	BOOST_AUTO_TEST_CASE(thrust_cpugpu_issue123_3D_complex) {
		using T         = thrust::complex<double>;
		auto const exts = multi::extensions_t<3>({1024, 1024, 100});

		std::cout << "| 3D `" << typeid(T).name() << "` max data size " << exts.num_elements() * sizeof(T) / 1073741824. << " GB | speed |\n|---|---|" << std::endl;

		multi::array<T, 3, test_allocator<T>> Devc(exts);
		multi::array<T, 3, test_allocator<T>> Dev2(exts);
		multi::array<T, 3>                    Host(exts);
		std::iota(Host.elements().begin(), Host.elements().end(), 12.);
		multi::array<T, 3> Hos2(exts);

		{
			Devc({0, 512}, {0, 512}, {0, 512}) = Host({0, 512}, {0, 512}, {0, 512});  // 0.002859s
		}
		{
			boost::timer::auto_cpu_timer t{""};
			Devc = Host;
			std::cout << "| contiguous host to devc | " << Host.num_elements() * sizeof(T) / (t.elapsed().wall / 1e9) / 1073741824. << " GB/sec |" << std::endl;
		}
		{
			boost::timer::auto_cpu_timer t{""};
			Devc.sliced(0, 512) = Host.sliced(0, 512);  //  0.005292s
			std::cout << "| sliced     host to devc | " << Host.sliced(0, 512).num_elements() * sizeof(T) / (t.elapsed().wall / 1e9) / 1073741824. << " GB/sec |" << std::endl;
		}
		{
			boost::timer::auto_cpu_timer t{""};
			Devc({0, 512}, {0, 512}, {0, 512}) = Host({0, 512}, {0, 512}, {0, 512});  // 0.002859s
			std::cout << "| strided    host to devc | " << Host({0, 512}, {0, 512}, {0, 512}).num_elements() * sizeof(T) / (t.elapsed().wall / 1e9) / 1073741824. << "GB/sec |" << std::endl;
		}
		{
			boost::timer::auto_cpu_timer t{""};
			Hos2 = Devc;
			std::cout << "| contiguous devc to host | " << Host.num_elements() * sizeof(T) / (t.elapsed().wall / 1e9) / 1073741824. << "GB/sec |" << std::endl;
			// BOOST_TEST(Hos2 == Host);
		}
		{
			boost::timer::auto_cpu_timer t{""};
			Hos2.sliced(0, 512) = Devc.sliced(0, 512);  //  0.005292s
			std::cout << "| sliced     devc to host | " << Host.sliced(0, 512).num_elements() * sizeof(T) / (t.elapsed().wall / 1e9) / 1073741824. << "GB/sec |" << std::endl;
			// BOOST_TEST(Hos2 == Host);
		}
		{
			boost::timer::auto_cpu_timer t{""};
			Hos2({0, 512}, {0, 512}, {0, 512}) = Devc({0, 512}, {0, 512}, {0, 512});  // 0.002859s
			std::cout << "| strided    devc to host | " << Host({0, 512}, {0, 512}, {0, 512}).num_elements() * sizeof(T) / (t.elapsed().wall / 1e9) / 1073741824. << "GB/sec |" << std::endl;
			// BOOST_TEST(Hos2 == Host);
		}
		{
			boost::timer::auto_cpu_timer t{""};
			Dev2 = Devc;
			cudaDeviceSynchronize();
			std::cout << "| contiguous devc to devc | " << Dev2.num_elements() * sizeof(T) / (t.elapsed().wall / 1e9) / 1073741824. << "GB/sec |" << std::endl;
			// BOOST_TEST(Dev2 == Devc);
		}
		{
			boost::timer::auto_cpu_timer t{""};
			auto                         Dev3 = Devc;
			cudaDeviceSynchronize();
			std::cout << "| copy_ctr   devc -> devc | " << Devc.num_elements() * sizeof(T) / (t.elapsed().wall / 1e9) / 1073741824. << "GB/sec |" << std::endl;
			// BOOST_TEST( Dev3 == Devc );
		}
		{
			boost::timer::auto_cpu_timer t{""};
			cudaMemcpy(raw_pointer_cast(Dev2.data_elements()), raw_pointer_cast(Devc.data_elements()), Devc.num_elements() * sizeof(T), cudaMemcpyDeviceToDevice);
			cudaDeviceSynchronize();
			std::cout << "| cudaMemcpy devc -> devc | " << Dev2.num_elements() * sizeof(T) / (t.elapsed().wall / 1e9) / 1073741824. << "GB/sec |" << std::endl;
			// BOOST_TEST( Dev2 == Devc );
		}
		{
			boost::timer::auto_cpu_timer t{""};
			cudaMemcpy(raw_pointer_cast(Dev2.data_elements()), raw_pointer_cast(Devc.data_elements()), Devc.num_elements() * sizeof(T), cudaMemcpyDeviceToDevice);
			cudaDeviceSynchronize();
			std::cout << "| cudaMemcpy devc -> devc | " << Dev2.num_elements() * sizeof(T) / (t.elapsed().wall / 1e9) / 1073741824. << "GB/sec |" << std::endl;
			// BOOST_TEST( Dev2 == Devc );
		}
		{
			boost::timer::auto_cpu_timer t{""};
			Dev2.sliced(0, 512) = Devc.sliced(0, 512);  //  0.005292s
			cudaDeviceSynchronize();
			std::cout << "| sliced     devc to devc | " << Dev2.sliced(0, 512).num_elements() * sizeof(T) / (t.elapsed().wall / 1e9) / 1073741824. << "GB/sec |" << std::endl;
			// BOOST_TEST(Dev2 == Devc);
		}
		{
			boost::timer::auto_cpu_timer t{""};
			Dev2({0, 512}, {0, 512}, {0, 512}) = Devc({0, 512}, {0, 512}, {0, 512});  // 0.002859s
			cudaDeviceSynchronize();
			std::cout << "| strided    devc to devc | " << Dev2({0, 512}, {0, 512}, {0, 512}).num_elements() * sizeof(T) / (t.elapsed().wall / 1e9) / 1073741824. << "GB/sec |" << std::endl;
			// BOOST_TEST(Dev2 == Devc);
		}
		{
			boost::timer::auto_cpu_timer t{""};
			Hos2 = Host;
			std::cout << "| contiguous host to host | " << Hos2.num_elements() * sizeof(T) / (t.elapsed().wall / 1e9) / 1073741824. << "GB/sec |" << std::endl;
		}
		{
			boost::timer::auto_cpu_timer t{""};
			Hos2.sliced(0, 512) = Host.sliced(0, 512);  //  0.005292s
			std::cout << "| sliced     host to host | " << Hos2.sliced(0, 512).num_elements() * sizeof(T) / (t.elapsed().wall / 1e9) / 1073741824. << "GB/sec |" << std::endl;
		}
		{
			boost::timer::auto_cpu_timer t{""};
			Hos2({0, 512}, {0, 512}, {0, 512}) = Host({0, 512}, {0, 512}, {0, 512});  // 0.002859s
			std::cout << "| strided    host to host | " << Hos2({0, 512}, {0, 512}, {0, 512}).num_elements() * sizeof(T) / (t.elapsed().wall / 1e9) / 1073741824. << "GB/sec |" << std::endl;
		}
		std::cout << "   " << std::endl;
	}
#endif

	return boost::report_errors();
}
#if 0

	#if 0
BOOST_AUTO_TEST_CASE_TEMPLATE(thrust_equality_1D_issue123, T, types_list) {
	multi::array<T, 1, test_allocator<T>> Devc(multi::extensions_t<1>{10240*10240});
	multi::array<T, 1, test_allocator<T>> Dev2(multi::extensions_t<1>{10240*10240});
	multi::array<T, 1>                    Host(multi::extensions_t<1>{10240*10240});
	std::iota(Host.elements().begin(), Host.elements().end(), 12.);
	multi::array<T, 1>                    Hos2(multi::extensions_t<1>{10240*10240});

	std::cout<<"| 1D `"<< typeid(T).name() <<"` total data size: "<< Host.num_elements()*sizeof(T) / 1073741824. <<" GB | speed |\n|---|---|\n";

	Devc = Host;
	Dev2 = Host;
	Hos2 = Host;
	{
		boost::timer::auto_cpu_timer t{""};
		BOOST_TEST( Devc == Dev2 );
		std::cout<<"| contiguous devc == devc | "<< Devc.num_elements()*sizeof(T) / (t.elapsed().wall/1e9) / 1073741824. << "GB/sec |\n";
	}
	{
		boost::timer::auto_cpu_timer t{""};
	//  BOOST_TEST( Devc.sliced(0, Devc.size()/2) == Dev2.sliced(0, Devc.size()/2) );
		BOOST_TEST( thrust::equal( Devc.sliced(0, Devc.size()/2).elements().begin(), Devc.sliced(0, Devc.size()/2).elements().end(), Dev2.sliced(0, Devc.size()/2).elements().begin() ) );
		std::cout<<"| sliced     devc == devc | "<< Devc.num_elements()*sizeof(T) / (t.elapsed().wall/1e9) / 1073741824. << "GB/sec |\n";
	}
	{
		boost::timer::auto_cpu_timer t{""};
	//  BOOST_TEST(Host == Hos2);
		BOOST_TEST( std::equal( Host.elements().begin(), Host.elements().end(), Hos2.elements().begin() ) );
		std::cout<<"| contiguous host == host | "<< Hos2.num_elements()*sizeof(T) / (t.elapsed().wall/1e9) / 1073741824. << "GB/sec |\n";
	}
	{
		boost::timer::auto_cpu_timer t{""};
	//  BOOST_TEST(Host.sliced(0, Devc.size()/2) == Hos2.sliced(0, Devc.size()/2) );
		BOOST_TEST( std::equal( Host.sliced(0, Host.size()/2).elements().begin(), Host.sliced(0, Host.size()/2).elements().end(), Hos2.sliced(0, Devc.size()/2).elements().begin() ) );
		std::cout<<"| sliced     host == host | "<< Hos2.num_elements()*sizeof(T) / (t.elapsed().wall/1e9) / 1073741824. << "GB/sec |\n";
	}

	std::cout<<"   "<<std::endl;
}


BOOST_AUTO_TEST_CASE(thrust_equality_2D_small_host_issue123) {
	multi::array<int, 2> Host = {{1, 2, 3}, {4, 5, 6}};
	multi::array<int, 2> Hos2 = {{1, 2, 3}, {4, 5, 6}};
	BOOST_TEST( Host.size() == 2 );

	BOOST_TEST( *Host().elements().begin()    == *Hos2().elements().begin()    );

	BOOST_TEST(  Host().elements().begin()[0] ==  Hos2().elements().begin()[0] );
	BOOST_TEST(  Host().elements().begin()[1] ==  Hos2().elements().begin()[1] );
	BOOST_TEST(  Host().elements().begin()[2] ==  Hos2().elements().begin()[2] );
	BOOST_TEST(  Host().elements().begin()[3] ==  Hos2().elements().begin()[3] );
	BOOST_TEST(  Host().elements().begin()[4] ==  Hos2().elements().begin()[4] );
	BOOST_TEST(  Host().elements().begin()[5] ==  Hos2().elements().begin()[5] );

	BOOST_TEST( *(Host().elements().end() - 1)    == *(Hos2().elements().end() - 1)   );
	BOOST_TEST( *(Host().elements().end() - 2)    == *(Hos2().elements().end() - 2)   );
	BOOST_TEST( *(Host().elements().end() - 3)    == *(Hos2().elements().end() - 3)   );

	BOOST_TEST(    std::equal(Host().elements().begin(), Host().elements().end(), Hos2().elements().begin()) );
	BOOST_TEST( thrust::equal(Host().elements().begin(), Host().elements().end(), Hos2().elements().begin()) );

//  BOOST_TEST( Host() == Hos2() );
}

BOOST_AUTO_TEST_CASE(thrust_equality_2D_small_gpu_issue123) {
	multi::array<int, 2> Host = {{1, 2, 3}, {4, 5, 6}};

	multi::array<int, 2, test_allocator<int>> Devc(Host.extensions()); Devc = Host;
	multi::array<int, 2, test_allocator<int>> Dev2(Host.extensions()); Dev2 = Host;
	BOOST_TEST( Dev2.size() == 2 );

	BOOST_TEST( thrust::equal(
		Devc().elements().begin(),
		Devc().elements().end()  , Dev2().elements().begin()
	));

	BOOST_TEST( thrust::equal(
		thrust::cuda::par,
		Devc().elements().begin(),
		Devc().elements().end()  , Dev2().elements().begin()
	));

	BOOST_TEST( thrust::equal(
		Devc.rotated().elements().begin(),
		Devc.rotated().elements().end()  , Dev2.rotated().elements().begin()
	));

	BOOST_TEST( thrust::equal(
		thrust::cuda::par,
		Devc.rotated().elements().begin(),
		Devc.rotated().elements().end()  , Dev2.rotated().elements().begin()
	));

	BOOST_TEST( multi::adl_equal(
		Devc.rotated().elements().begin(),
		Devc.rotated().elements().end()  , Dev2.rotated().elements().begin()
	));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(thrust_equality_2D_issue123, T, types_list) {
	multi::extensions_t<2> x({10240, 10240});
	multi::array<T, 2, test_allocator<T>> Devc(x);
	multi::array<T, 2, test_allocator<T>> Dev2(x);
	multi::array<T, 2>                    Host(x); std::iota(Host.elements().begin(), Host.elements().end(), 12.);
	multi::array<T, 2>                    Hos2(x);

	std::cout<<"| 2D `"<< typeid(T).name() <<"` max data size "<< Host.num_elements()*sizeof(T) / 1073741824. <<" GB | speed |\n|---|---|\n";

	Devc = Host;
	Dev2 = Host;
	Hos2 = Host;
	{
		boost::timer::auto_cpu_timer t{""};
		BOOST_TEST(Host == Hos2);
		std::cout<<"| contiguous host == host | "<< Hos2.num_elements()*sizeof(T) / (t.elapsed().wall/1e9) / 1073741824. << "GB/sec |\n";
	}
	{
		boost::timer::auto_cpu_timer t{""};
		BOOST_TEST(Host.sliced(0, Host.size()/2) == Hos2.sliced(0, Host.size()/2));
	//  BOOST_TEST( std::equal(Host.sliced(0, 5120).elements().begin(), Host.sliced(0, 5120).elements().end(), Hos2.sliced(0, 5120).elements().begin()) );
		std::cout<<"| sliced     host == host | "<< Hos2.num_elements()*sizeof(T) / (t.elapsed().wall/1e9) / 1073741824. << "GB/sec |\n";
	}
	{
		boost::timer::auto_cpu_timer t{""};
		BOOST_TEST(Host({0, Host.size()/2},{0, Host.size()/2}) == Hos2({0, Hos2.size()/2},{0, Hos2.size()/2}));
	//  BOOST_TEST( std::equal(Host({0, 5120},{0, 5120}).elements().begin(), Host({0, 5120},{0, 5120}).elements().end(), Hos2({0, 5120},{0, 5120}).elements().begin()) );
		std::cout<<"| strided    host == host | "<< Hos2({0, Host.size()/2},{0, Host.size()/2}).num_elements()*sizeof(T) / (t.elapsed().wall/1e9) / 1073741824. << "GB/sec |\n";
	}
	{
		boost::timer::auto_cpu_timer t{""};
		BOOST_TEST(Devc == Dev2);
		std::cout<<"| contiguous devc == devc | "<< Devc.num_elements()*sizeof(T) / (t.elapsed().wall/1e9) / 1073741824. << "GB/sec |\n";
	}
	{
		boost::timer::auto_cpu_timer t{""};
		BOOST_TEST(Devc.sliced(0, Devc.size()/2) == Dev2.sliced(0, Dev2.size()/2));
		std::cout<<"| sliced     devc == devc | "<< Devc.sliced(0, Devc.size()/2).num_elements()*sizeof(T) / (t.elapsed().wall/1e9) / 1073741824. << "GB/sec |\n";
	}
	{
		boost::timer::auto_cpu_timer t{""};
		BOOST_TEST(Devc({0, Devc.size()/2},{0, Devc.size()/2}) == Dev2({0, Dev2.size()/2},{0, Dev2.size()/2}));
		std::cout<<"| strided    devc == devc | "<< Devc({0, Devc.size()/2},{0, Devc.size()/2}).num_elements()*sizeof(T) / (t.elapsed().wall/1e9) / 1073741824. << "GB/sec |\n";
	}
	std::cout<<"   "<<std::endl;
}

		#if 1
BOOST_AUTO_TEST_CASE_TEMPLATE(thrust_equality_issue123_3D, T, types_list) {
	multi::array<T, 3, test_allocator<T>> Devc({1024, 1024, 100});
	multi::array<T, 3, test_allocator<T>> Dev2({1024, 1024, 100});
	multi::array<T, 3>                    Host({1024, 1024, 100}); std::iota(Host.elements().begin(), Host.elements().end(), 12.);
	multi::array<T, 3>                    Hos2({1024, 1024, 100});

	std::cout<<"| 3D `"<< typeid(T).name() <<"` max data size "<< Host.num_elements()*sizeof(T) / 1073741824. <<" GB | speed |\n|---|---|\n";

	Devc = Host;
	Dev2 = Host;
	Hos2 = Host;
	{
		boost::timer::auto_cpu_timer t{""};
		BOOST_TEST( Devc == Dev2 );
		std::cout<<"| contiguous devc == devc | "<< Devc.num_elements()*sizeof(T) / (t.elapsed().wall/1e9) / 1073741824. << " GB/sec |\n";
	}
	{
		boost::timer::auto_cpu_timer t{""};
		BOOST_TEST( Devc.sliced(0, 512) == Dev2.sliced(0, 512) );
		std::cout<<"| sliced     devc == devc | "<< Dev2.sliced(0, 512).num_elements()*sizeof(T) / (t.elapsed().wall/1e9) / 1073741824. << " GB/sec |\n";
	}
	{
		boost::timer::auto_cpu_timer t{""};
		BOOST_TEST( Devc({0, 512}, {0, 512}, {0, 512}) == Dev2({0, 512}, {0, 512}, {0, 512}) );
		std::cout<<"| strided    devc == devc | "<< Dev2({0, 512},{0, 512}, {0, 512}).num_elements()*sizeof(T) / (t.elapsed().wall/1e9) / 1073741824. << "GB/sec |\n";
	}
	{
		boost::timer::auto_cpu_timer t{""};
		BOOST_TEST( Host == Hos2 );
		std::cout<<"| contiguous host == host | "<< Hos2.num_elements()*sizeof(T) / (t.elapsed().wall/1e9) / 1073741824. << " GB/sec |\n";
	}
	{
		boost::timer::auto_cpu_timer t{""};
	//  BOOST_TEST( Host.sliced(0, 512) == Hos2.sliced(0, 512) );
		BOOST_TEST( std::equal( Host.sliced(0, 512).elements().begin(), Host.sliced(0, 512).elements().end(),  Hos2.sliced(0, 512).elements().begin() ) );
		std::cout<<"| sliced     host == host | "<< Hos2.sliced(0, 512).num_elements()*sizeof(T) / (t.elapsed().wall/1e9) / 1073741824. << " GB/sec |\n";
	}
	{
		boost::timer::auto_cpu_timer t{""};
	//  BOOST_TEST( Host({0, 512}, {0, 512}, {0, 512}) == Hos2({0, 512}, {0, 512}, {0, 512}) );
		BOOST_TEST( std::equal( Host({0, 512}, {0, 512}, {0, 512}).elements().begin(), Host({0, 512}, {0, 512}, {0, 512}).elements().end(), Hos2({0, 512}, {0, 512}, {0, 512}).elements().begin() ) );
		std::cout<<"| strided    host == host | "<< Hos2({0, 512},{0, 512}, {0, 512}).num_elements()*sizeof(T) / (t.elapsed().wall/1e9) / 1073741824. << "GB/sec |\n";
	}
	std::cout<<"   "<<std::endl;
}
		#endif
	#endif

	#if 0
namespace inq {
	using complex = thrust::complex<double>;
}

BOOST_AUTO_TEST_CASE(thrust_complex_cached_1D) {
	using T = inq::complex;
	multi::array<T, 1, boost::multi::memory::cuda::cached::allocator<T> > aa(10, T{1., 1.});
	multi::array<T, 1, boost::multi::memory::cuda::cached::allocator<T> > bb(10, T{2., 2.});

	bb = aa;

	BOOST_TEST(( bb[0] == T{1., 1.} ));
}

BOOST_AUTO_TEST_CASE(thrust_complex_cached_without_values_1D) {
	using T = inq::complex;
	multi::array<T, 1, boost::multi::memory::cuda::cached::allocator<T> > aa(10);
	multi::array<T, 1, boost::multi::memory::cuda::cached::allocator<T> > bb(10);
	BOOST_TEST( aa.size() == 10 );
	BOOST_TEST( bb.size() == 10 );

	bb = aa;

	BOOST_TEST(( bb[0] == aa[0] ));
}

BOOST_AUTO_TEST_CASE(thrust_complex_cached_2D) {
	using T = inq::complex;
	multi::array<T, 2, boost::multi::memory::cuda::cached::allocator<T> > aa({10, 20}, T{1., 1.});
	multi::array<T, 2, boost::multi::memory::cuda::cached::allocator<T> > bb({10, 20}, T{2., 2.});

	bb = aa;

	BOOST_TEST(( bb[0][0] == T{1., 1.} ));
}

BOOST_AUTO_TEST_CASE(thrust_complex_cached_without_values_2D) {
	using T = inq::complex;
	multi::array<T, 2, boost::multi::memory::cuda::cached::allocator<T> > aa({10, 20});
	multi::array<T, 2, boost::multi::memory::cuda::cached::allocator<T> > bb({10, 20});
	BOOST_TEST( aa.size() == 10 );
	BOOST_TEST( bb.size() == 10 );

	bb = aa;

	BOOST_TEST(( bb[0][0] == aa[0][0] ));
}

BOOST_AUTO_TEST_CASE(array) {

//{
//  multi::thrust::cuda::array<double, 2> C({2, 3});

//  C[0][0] = 0. ;
//  C[1][1] = 11.;
//  BOOST_TEST_REQUIRE( C[1][1] == 11. );
//}

//{
//  multi::array<double, 2> const H = {
//      {00., 01., 02.},
//      {10., 11., 12.},
//  };

//  BOOST_TEST_REQUIRE( H[1][1] == 11. );

//  {
//      multi::thrust::cuda::array<double, 2> C(H.extensions());
//      BOOST_TEST( C.num_elements() == H.num_elements() );

//      thrust::copy_n(H.data_elements(), H.num_elements(), C.data_elements());
//      BOOST_TEST_REQUIRE( C[1][1] == 11. );
//      BOOST_TEST( C == H );
//  }
//  {
//      multi::thrust::cuda::array<double, 2> C(H.extensions());
//      BOOST_TEST( C.num_elements() == H.num_elements() );

//      std::copy_n(H.data_elements(), H.num_elements(), C.data_elements());
//      BOOST_TEST_REQUIRE( C[1][1] == 11. );
//      BOOST_TEST( C == H );
//  }
//  {
//      multi::thrust::cuda::array<double, 2> C(H.extensions());
//      BOOST_TEST( C.num_elements() == H.num_elements() );

//      std::uninitialized_copy_n(H.data_elements(), H.num_elements(), C.data_elements());
//      BOOST_TEST_REQUIRE( C[1][1] == 11. );
//      BOOST_TEST( C == H );
//  }
//  {
//      multi::thrust::cuda::array<double, 2> C(H.extensions());
//      BOOST_TEST( C.num_elements() == H.num_elements() );

//      thrust::uninitialized_copy_n(H.data_elements(), H.num_elements(), C.data_elements());
//      BOOST_TEST_REQUIRE( C[1][1] == 11. );
//      BOOST_TEST( C == H );
//  }
//  {
//      multi::thrust::cuda::array<double, 2> C(H.extensions());
//      BOOST_TEST( C.extensions() == H.extensions() );
//      thrust::copy_n(H.begin(), H.size(), C.begin());
//      BOOST_TEST( C == H );
//  }
//  {
//      multi::thrust::cuda::array<double, 2> C(H.extensions());
//      BOOST_TEST( C.extensions() == H.extensions() );
//      std::copy_n(H.begin(), H.size(), C.begin());
//      BOOST_TEST( C == H );
//  }
//  {
//      multi::thrust::cuda::array<double, 2> C(H.extensions());
//      C = H;
//      BOOST_TEST( C == H );
//  }
//  {
//      multi::thrust::cuda::array<double, 2> C = H;
//      BOOST_TEST( C == H );
//  }
//}

}

return boost::report_errors();}

	#endif
#endif