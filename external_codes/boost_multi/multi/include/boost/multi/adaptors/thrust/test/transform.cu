// Copyright 2026 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/multi/adaptors/thrust.hpp>
#include <boost/multi/array.hpp>

#include <thrust/complex.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/random.h>
#include <thrust/tabulate.h>

#include <boost/core/lightweight_test.hpp>

namespace multi = boost::multi;

template<typename T>
T* start_lifetime_as_array_approx(void* p, std::size_t n) {
	// memmove to self is a "blessed" operation that
	// implicitly starts lifetimes for implicit-lifetime types.
	std::memmove(p, p, sizeof(T) * n);
	return reinterpret_cast<T*>(p);
}

auto main() -> int {
	namespace multi = boost::multi;

	multi::thrust::device_array<thrust::complex<double>, 3> olap({5, 5, 3}, thrust::complex<double>{1.0, 2.0});

	// ::thrust::tabulate(olap.elements().begin(), olap.elements().end(), [seed = 42] __device__ (multi::index i) {
	// 	thrust::default_random_engine rng(seed);
	// 	rng.discard(i);
	// 	thrust::uniform_real_distribution<double> dist(-1.0, 1.0);
	// 	return thrust::complex<double>(dist(rng), dist(rng));
	// });

	multi::thrust::device_array<double, 2> vel_gold;

	{
		multi::array<thrust::complex<double>, 3> olap_cpu(olap);
		multi::array<double, 2>                  vel({5, 5});

		std::transform(
			olap_cpu.flatted().begin(), olap_cpu.flatted().end(),
			vel.flatted().begin(),
			[](auto const& e) { return norm(e[0]) + norm(e[1]) + norm(e[2]); }
		);

		vel_gold = vel;

		BOOST_TEST( std::abs(vel[2][3] - vel_gold[2][3]) < 1e-12  );
	}
	// {
	// 	multi::array<double, 2> vel({5, 5});

	// 	thrust::transform(
	// 		olap.reinterpret_array_cast<std::array<thrust::complex<double>, 3>>().flatted().begin(),
	// 		olap.reinterpret_array_cast<std::array<thrust::complex<double>, 3>>().flatted().end(),
	// 		vel.flatted().begin(),
	// 		[] __device__ (auto const& e) { return norm(e[0]) + norm(e[1]) + norm(e[2]); }
	// 	);

	// 	vel_gold = vel;

	// 	BOOST_TEST( std::abs(vel[2][3] - vel_gold[2][3]) < 1e-12  );
	// }
	{
		multi::thrust::universal_array<double, 2> vel({5, 5});

		thrust::transform(
			vel.elements().extent().begin(),
			vel.elements().extent().end(),
			vel.elements().begin(),
			[olap_base = olap.flatted().home()] __device__(int mm) {
				return thrust::norm(static_cast<thrust::complex<double>>(olap_base[mm][0])) + thrust::norm(static_cast<thrust::complex<double>>(olap_base[mm][1])) + thrust::norm(static_cast<thrust::complex<double>>(olap_base[mm][2]));
			}
		);

		BOOST_TEST( std::abs(vel[2][3] - vel_gold[2][3]) < 1e-12  );
	}
#ifndef _WIN32
	{
		multi::thrust::device_array<double, 2> vel({5, 5});

		thrust::transform(
			thrust::make_counting_iterator(olap.flatted().begin()),
			thrust::make_counting_iterator(olap.flatted().end()),
			vel.elements().begin(),
			[] __device__(auto e) {
				auto const& ev = *e;
				return norm(static_cast<thrust::complex<double>>(ev[0])) + norm(static_cast<thrust::complex<double>>(ev[1])) + norm(static_cast<thrust::complex<double>>(ev[2]));
			}
		);

		BOOST_TEST( std::abs(vel[2][3] - vel_gold[2][3]) < 1e-12  );
	}
#endif
	{
		multi::thrust::device_array<int, 1> A({100}, 3);
		multi::thrust::device_array<int, 1> B({100}, 2);
		multi::thrust::device_array<int, 1> C({100}, -1);

		thrust::transform(
			thrust::make_zip_iterator(A.begin(), B.begin()),
			thrust::make_zip_iterator(A.end(), B.end()),
			C.begin(),
			[] __device__(auto const& ab) {
				using std::get;  // can't use structured bindings with thrust::tuple
				return static_cast<int>(get<0>(ab)) + static_cast<int>(get<1>(ab));
			}
		);

		BOOST_TEST( C[21] == 5 );

		multi::array<int, 1> C_cpu(C);
		std::cout << "C_cpu[21] = " << C_cpu[21] << '\n';
		BOOST_TEST( C_cpu[21] == 5 );
	}
	{
		multi::thrust::device_array<int, 2> AB({100, 2}, -1);

		auto&& A = AB(multi::_, 0);
		auto&& B = AB(multi::_, 1);

		thrust::fill(A.begin(), A.end(), 3);
		thrust::fill(B.begin(), B.end(), 2);

		{
			multi::thrust::device_array<int, 1> C({100}, -1);

			thrust::transform(
				thrust::make_zip_iterator(A.begin(), B.begin()),
				thrust::make_zip_iterator(A.end(), B.end()),
				C.begin(),
				[] __device__(auto const& ab) {
					using std::get;  // can't use structured bindings with thrust::tuple
					return get<0>(ab) + get<1>(ab);
				}
			);

			BOOST_TEST( C[21] == 5 );

			multi::array<int, 1> C_cpu(C);
			std::cout << "C_cpu[21] = " << C_cpu[21] << '\n';
			BOOST_TEST( C_cpu[21] == 5 );
		}
		{
			using iter1d = decltype(std::declval<multi::thrust::device_array<int, 1>>().begin());
			static_assert(std::is_same_v<typename thrust::iterator_system<iter1d>::type, thrust::device_system_tag>, "1D iterator should be device");

			using iter2d = decltype(std::declval<multi::thrust::device_array<int, 2>>().begin());
			static_assert(std::is_same_v<typename thrust::iterator_system<iter2d>::type, thrust::device_system_tag>, "2D iterator should be device");

#ifndef _WIN32
			multi::thrust::device_array<int, 1> C({100}, -1);

			thrust::transform(
				AB.begin(),
				AB.end(),
				C.begin(),
				[] __device__(auto const& ab) {
					return ab[0] + ab[1];
				}
			);

			BOOST_TEST( C[21] == 5 );

			multi::array<int, 1> C_cpu(C);
			BOOST_TEST( C_cpu[21] == 5 );
#endif
		}
		{
#ifndef _WIN32
			multi::thrust::device_array<int, 1> C({100}, -1);

			thrust::transform(
				thrust::cuda::par,
				AB.begin(),
				AB.end(),
				C.begin(),
				[] __device__(auto const& ab) {
					return ab[0] + ab[1];
				}
			);

			BOOST_TEST( C[21] == 5 );

			multi::array<int, 1> C_cpu(C);
			BOOST_TEST( C_cpu[21] == 5 );
#endif
		}
		{
#ifndef _WIN32
			multi::thrust::device_array<int, 1> C({100}, -1);

			thrust::transform(
				thrust::cuda::par,
				AB.begin(),
				AB.end(),
				C.begin(),
				[] __device__ (auto const& ab) {
					return ab[0] + ab[1];
				}
			);

			BOOST_TEST( C[21] == 5 );

			multi::array<int, 1> C_cpu(C);
			BOOST_TEST( C_cpu[21] == 5 );
#endif
		}

	}

	return boost::report_errors();
}
