// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2020-2023 Alfredo A. Correa

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi FFTW adaptor"

#include<boost/test/unit_test.hpp>

#include <multi/adaptors/fftw.hpp>

#include<chrono>
#include<complex>
#include<iostream>
#include<random>

namespace multi = boost::multi;

template<>
inline constexpr bool multi::force_element_trivial_default_construction<std::complex<double>> = true;

namespace utf = boost::unit_test::framework;

using namespace std::string_literals;  // NOLINT(build/namespaces) for ""s

using fftw_fixture = multi::fftw::environment;
BOOST_TEST_GLOBAL_FIXTURE( fftw_fixture );

class watch : private std::chrono::high_resolution_clock {
	std::string label;
	time_point start = now();

 public:
	template<class String>
	explicit watch(String&& label) : label{std::forward<String>(label)} {}  // NOLINT(fuchsia-default-arguments-calls)

	watch(watch const&) = delete;
	watch(watch&&) = delete;

	auto operator=(watch const&) = delete;
	auto operator=(watch&&) = delete;

	auto elapsed_sec() const {return std::chrono::duration<double>(now() - start).count();}
	~watch() { std::cerr<< label <<": "<< elapsed_sec() <<" sec"<<std::endl; }
};

template<class T, multi::dimensionality_type D> using marray = multi::array<T, D>;
constexpr auto exts = multi::extensions_t<4>({6, 12, 24, 12});

BOOST_AUTO_TEST_CASE(fft_combinations, *boost::unit_test::tolerance(0.00001) ) {
	using complex = std::complex<double>;

	auto const in = [] {
		marray<complex, 4> ret(exts);
		std::generate(ret.data_elements(), ret.data_elements() + ret.num_elements(),
			[eng = std::default_random_engine {std::random_device {}()},
				uniform_01 = std::uniform_real_distribution<>{}]() mutable {
				return complex{uniform_01(eng), uniform_01(eng)};
			});
		return ret;
	}();

	std::vector<std::array<bool, 4>> const which_cases = {  // std::vector NOLINT(fuchsia-default-arguments-calls)
		{false, true , true , true },
		{false, true , true , false},
		{true , false, false, false},
		{true , true , false, false},
		{false, false, true , false},
		{false, false, false, false},
	};

	using std::cout;
	for(auto which : which_cases) {  // NOLINT(altera-unroll-loops)
		cout<<"case ";
		copy(begin(which), end(which), std::ostream_iterator<bool>{cout, ", "});
		cout<<"\n";

		marray<complex, 4> out = in;
		// {
		//  watch const unnamed{"cpu_oplac %ws wall, CPU (%p%)\n"s};
		//  multi::fftw::dft_forward(which, in, out);
		// }
		{
			auto const pln = multi::fftw::plan::forward(which, in.base(), in.layout(), out.base(), out.layout());
			watch const unnamed{"cpu_oplac planned %ws wall, CPU (%p%)\n"s};
			pln.execute(in.base(), out.base());
		}
		{
			auto in_rw = in;
			watch const unnamed{"cpu_iplac %ws wall, CPU (%p%)\n"s};
			multi::fftw::dft_forward(which, in_rw, in_rw);
		}
		{
			auto in_rw = in;
			auto const pln = multi::fftw::plan::forward(which, in_rw.base(), in_rw.layout(), in_rw.base(), in_rw.layout());
			watch const unnamed{"cpu_iplac planned %ws wall, CPU (%p%)\n"s};
			pln.execute(in_rw.base(), in_rw.base());
		}
		{
			auto in_rw = in;
			auto const pln = multi::fftw::plan::forward(which, in_rw.base(), in_rw.layout(), in_rw.base(), in_rw.layout());
			watch const unnamed{"cpu_iplac planned measured %ws wall, CPU (%p%)\n"s};
			pln.execute(in_rw.base(), in_rw.base());
		}
		// {
		//  watch const unnamed{"cpu_alloc %ws wall, CPU (%p%)\n"s};
		//  auto out_cpy = multi::fftw::dft_forward(which, in);
		//  BOOST_TEST(abs(out_cpy[5][4][3][1] - out[5][4][3][1]) == 0.);
		// }
		// {
		//  auto in_rw = in;
		//  watch const unnamed{"cpu_move %ws wall, CPU (%p%)\n"s};
		//  auto out_cpy = multi::fftw::dft_forward(which, std::move(in_rw));
		//  BOOST_TEST(abs(out_cpy[5][4][3][1] - out[5][4][3][1]) == 0.);
		// }
	}
}

BOOST_AUTO_TEST_CASE(fftw_4D_power_benchmark, *boost::unit_test::enabled() ) {
	using complex = std::complex<double>;
	namespace fftw = multi::fftw;

	marray<complex, 4> in(exts);
	std::iota(in.data_elements(), in.data_elements() + in.num_elements(), 1.2);

	BOOST_REQUIRE(in[0][0][0][0] == 1.2);
	std::array<bool, 4> which = {false, true, true, true};
	[&, unnamed = watch{utf::current_test_case().full_name()+" inplace FTTT"s}] {
		fftw::dft(which, in, fftw::forward);
	}();
	[&, unnamed = watch{utf::current_test_case().full_name()+" inplace FTTT"s}] {
		fftw::dft(which, in, fftw::forward);
	}();
	auto in0000 = in[0][0][0][0];
	BOOST_REQUIRE(in0000 != 1.2);

	marray<complex, 4> out(exts);
	[&, unnamed = watch{utf::current_test_case().full_name()+" outofplace FTTT"s}] {
		fftw::dft(which, in, out, fftw::forward);
	}();
	[&, unnamed = watch{utf::current_test_case().full_name()+" outofplace FTTT"s}] {
		fftw::dft(which, in, out, fftw::forward);
	}();
	[&, unnamed = watch{utf::current_test_case().full_name()+" outofplace FTTT"s}] {
		fftw::dft(which, in, out, fftw::forward);
	}();
	[&, unnamed = watch{utf::current_test_case().full_name()+" outofplace+alloc FTTT"s}] {
		marray<complex, 4> out2(exts);
		fftw::dft(which, in, out2, fftw::forward);
	}();
	[&, unnamed = watch{utf::current_test_case().full_name()+" outofplace+alloc FTTT"s}] {
		marray<complex, 4> out2(exts);
		fftw::dft(which, in, out2, fftw::forward);
	}();
	BOOST_REQUIRE(in0000 == in[0][0][0][0]);
}

BOOST_AUTO_TEST_CASE(fftw_4D_power_benchmark_syntax) {
	std::vector<std::array<bool, 4>> const which_cases = {  // std::vector NOLINT(fuchsia-default-arguments-calls)
		{false, true , true , true },
		{false, true , true , false},
		{true , false, false, false},
		{true , true , false, false},
		{false, false, true , false},
		{false, false, false, false},
	};
	using complex = std::complex<double>;

	auto const in = [] {
		marray<complex, 4> ret(exts);
		std::generate(ret.data_elements(), ret.data_elements() + ret.num_elements(),
			[eng = std::default_random_engine {std::random_device {}()},
				uniform_01 = std::uniform_real_distribution<>{}]() mutable{
				return complex{uniform_01(eng), uniform_01(eng)};
			});
		return ret;
	}();

	auto io = in; (void)io;
	BOOST_REQUIRE( io.extensions() == in.extensions() );

	// namespace fftw = multi::fftw;
	// using clock = std::chrono::high_resolution_clock;
	// {
	//  auto const tick = clock::now();
	//  marray<complex, 4> out(exts);
	//  out = multi::fftw::ref(in)(fftw::none, fftw::forward, fftw::forward, fftw::forward);
	//  BOOST_REQUIRE( out.extensions() == in.extensions() );
	//  auto time = std::chrono::duration<double>(clock::now() - tick);
	//  std::cout<<"allocate and copy assign (out-of-place fft) : "<< time.count() <<std::endl;
	// }
	// {
	//  auto const tick = clock::now();
	//  auto const out = +multi::fftw::ref(in)(fftw::none, fftw::forward, fftw::forward, fftw::forward);
	//  BOOST_REQUIRE( out.extensions() == in.extensions() );
	//  auto time = std::chrono::duration<double>(clock::now() - tick);
	//  std::cout<<"copy construct (out-of-place fft) : "<< time.count() <<std::endl;
	// }
	// {
	//  auto const tick = clock::now();
	//  io = multi::fftw::ref(io)(fftw::none, fftw::forward, fftw::forward, fftw::forward);
	//  BOOST_REQUIRE( io.extensions() == in.extensions() );
	//  auto time = std::chrono::duration<double>(clock::now() - tick);
	//  std::cout<<"self copy assign (in-place fft) : "<< time.count() <<std::endl;
	// }
	// {
	//  auto const tick = clock::now();
	//  marray<complex, 4> const out = multi::fftw::move(io)(fftw::none, fftw::forward, fftw::forward, fftw::forward);
	//  BOOST_REQUIRE( io.is_empty() );
	//  auto time = std::chrono::duration<double>(clock::now() - tick);
	//  std::cout<<"move construct (in-place fft) : "<< time.count() <<std::endl;
	// }
}
