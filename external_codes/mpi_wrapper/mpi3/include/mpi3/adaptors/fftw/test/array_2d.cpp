// Copyright 2019-2024 Alfredo A. Correa

#include <mpi3/adaptors/fftw.hpp>
#include <mpi3/main_environment.hpp>
#include <mpi3/ostream.hpp>

#include <complex>  // for std::norm

template<class M> auto power(M const& elem) -> decltype(std::norm(elem)) { return std::norm(elem); }

template<class M, class = std::enable_if_t<(M::rank::value >= 1)>>  // DELETE((M::rank::value < 1))>
auto power(M const& array) {
	return accumulate(begin(array), end(array), 0.0, [](auto const& alpha, auto const& omega) { return alpha + power(omega); });
}

struct sum_power {
	template<class A, class B> auto operator()(A const& alpha, B const& omega) const { return alpha + power(omega); }
};

// template<class Array>
// void chop(Array&& arr) {
//  std::replace_if(arr.elements().begin(), arr.elements().end(), [](auto const& e) {std::fabs(e) < 1.0e-30}, 0.0);
//  // for(auto& e : arr.elements()) {
//  //  if(std::fabs(e) < 1.0e-30) {
//  //    e = 0.0;
//  //  }
//  // }
// }

template<class Array>
void mpi_fill(Array&& arr) {
	auto [is, js] = arr.local_cutout().extensions();
	std::for_each(is.begin(), is.end(), [&, js = js](auto i) {
		std::for_each(js.begin(), js.end(), [&](auto j) {
			arr.local_cutout()[i][j] = std::complex<double>{static_cast<double>(i + j), static_cast<double>(i + 2 * j)};
		});
	});
}

template<class Array>
void mpi_print(Array const& arr, boost::mpi3::communicator& comm, std::string const& /*msg*/ = "") {
	boost::mpi3::ostream ccout{comm, std::cout};

	ccout << "rank=" << comm.rank() << " count=" << arr.local_count() << '\n';
	auto [is, js] = arr.local_cutout().extensions();
	std::for_each(is.begin(), is.end(), [&, js = js](auto i) {
		std::for_each(js.begin(), js.end(), [&](auto j) {
			ccout << arr.local_cutout()[i][j] << " ";
		});
		ccout << '\n';
	});
}

namespace mpi3 = boost::mpi3;

auto mpi3::main(int /*argc*/, char** /*argv*/, boost::mpi3::environment& env) -> int try {
	auto world = env.world();

	boost::mpi3::fftw::environment fftwenv;

	boost::mpi3::fftw::array<std::complex<double>, 2, boost::mpi3::fftw::local_2d> G({6, 6}, 0.0, world);
	mpi_fill(G);

	boost::mpi3::fftw::array<std::complex<double>, 2, boost::mpi3::fftw::local_2d> F({6, 6}, 0.0, world);
	dft_forward(G, F);

	boost::multi::array<std::complex<double>, 2> const g{G};
	boost::multi::array<std::complex<double>, 2> f(g.extensions());

	boost::multi::fftw::dft_forward({true, true}, g, f);

	boost::multi::array<std::complex<double>, 2> const ff{F};
	assert( ff == f );


	// boost::mpi3::fftw::array<std::complex<double>, 2, boost::mpi3::fftw::local_2d_many> G_many({6, 6}, 0.0, world);

	// mpi_print(G_many, world);

	// dft(G, G_many);

	// mpi_print(G_many, world);

//  G_many = G;

	if(world.rank() == 0) {
		// assert( G_many.local_cutout()[2][2] == std::complex<double>(4.0, 6.0) );
	}

	// multi::array<std::complex<double>, 2> g{G};

	// if(world.rank() == 0) {
	//  std::cout << "gathered power " << power(g) << std::endl;
	//  auto [is, js] = g.extensions();
	//  for(auto i : is) {
	//    for(auto j : js) {
	//      std::cout << g[i][j] << ",";
	//    }
	//    std::cout << std::endl;
	//  }
	// }
	world.barrier();

	// boost::mpi3::fftw::array<std::complex<double>, 2> F({6, 6}, 0.0, world);
	// F.scatter(g);

	// multi::array<std::complex<double>, 2> f{F};
	// assert(g == f);

	// auto F2 = boost::mpi3::fftw::array<std::complex<double>, 2>::from_scatter(g);
//  auto F3 = boost::mpi3::scatter(g);

//  multi::array<std::complex<double>, 2> f3{F3};

	// assert(f3 == g);

	// multi::array<std::complex<double>, 2> g_transformed(g.extensions());
	// boost::multi::fftw::dft_forward({true, true}, g, g_transformed);

	// if(world.rank() == 0) {
	//  std::cout << "g_transformed power " << power(g_transformed) / g_transformed.num_elements() << std::endl;
	//  auto [is, js] = g_transformed.extensions();
	//  for(auto i : is) {
	//    for(auto j : js) {
	//      std::cout << g_transformed[i][j] << ",";
	//    }
	//    std::cout << std::endl;
	//  }
	// }
	// world.barrier();

	// multi::fftw::mpi::array<std::complex<double>, 2> G_transformed({6, 6}, 0.0, &world);
	// boost::multi::fftw::mpi::dft_forward(G, G_transformed);
	// chop(G_transformed);

	// boost::multi::fftw::mpi::dft_forward(G, G);
	// chop(G);

	// mpi_print(G, world, "G_transformed");

	// multi::array<std::complex<double>, 2> g_mpi_transformed{G};
	// chop(g_mpi_transformed);

	// if(world.rank() == 0) {
	//  std::cout << "g_mpi_transformed power " << power(g_mpi_transformed) / g_mpi_transformed.num_elements() << std::endl;
	//  auto [is, js] = g_mpi_transformed.extensions();
	//  for(auto i : is) {
	//    for(auto j : js) {
	//      std::cout << g_mpi_transformed[i][j] << ",";
	//    }
	//    std::cout << std::endl;
	//  }
	// }
	// world.barrier();

	// assert(g_mpi_transformed == g_transformed);

	// world.barrier();
	// mpi_print(G, world, "G transformed already");

	// world.barrier();
	// mpi_print(G_transformed, world, "G_transformed");

	return 0;
} catch(...) {
	return 1;
}
