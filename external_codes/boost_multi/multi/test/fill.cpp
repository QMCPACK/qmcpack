// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// © Alfredo A. Correa 2019-2021

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi fill"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include "../array.hpp"

#include<limits>
#include<random>
#include<type_traits> // enable_if_t

// from Howard Hinnart hash
auto fnv1a(void const* key, std::size_t len, std::size_t h) noexcept {  // NOLINT(bugprone-easily-swappable-parameters)
	auto const *p = static_cast<unsigned char const*>(key);
	unsigned char const* const e = p + len; // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic): low level
	for(; p < e; ++p) {  // NOLINT(altera-id-dependent-backward-branch,cppcoreguidelines-pro-bounds-pointer-arithmetic): low level
		h = (h ^ *p) * 1099511628211U; // prime
	}
	return h;
}

auto fnv1a(void const* key, std::size_t len) noexcept;
auto fnv1a(void const* key, std::size_t len) noexcept {
	return fnv1a(key, len, 14695981039346656037U);
}

class fnv1a_t{
	std::size_t h = 14695981039346656037U; // offset

 public:
	using result_type = std::size_t;
	static constexpr auto min() {return std::numeric_limits<result_type>::min();}
	static constexpr auto max() {return std::numeric_limits<result_type>::max();}
	void operator()(void const* key, std::size_t len) noexcept {h = fnv1a(key, len, h);}
	template<class T, std::enable_if_t<std::is_fundamental<T>{}, int> = 0>
	auto operator()(T const& t) noexcept -> decltype(auto) {operator()(&t, sizeof(t)); return *this;}
//	result_type operator()() && noexcept{return h;}
	auto operator()() const& noexcept {return h;}
//	explicit operator result_type() && noexcept{return h;}
	explicit operator result_type() const& noexcept {return h;}
};

std::random_device r;

BOOST_AUTO_TEST_CASE(fill_1d) {
	namespace multi = boost::multi;
	{
		multi::array<double, 1> d1D(multi::extensions_t<1>{multi::iextension{10}});
		static_assert( std::is_same<std::iterator_traits<decltype(begin(d1D))>::value_type, double>{}, "!");

		using std::copy;
		copy(begin(extension(d1D)), end(extension(d1D)), begin(d1D));
		BOOST_REQUIRE( d1D[0] == 0. );
		BOOST_REQUIRE( d1D[1] == 1. );
		BOOST_REQUIRE( d1D[9] == 9. );

		d1D.assign(extension(d1D));
		BOOST_REQUIRE( d1D[0] == 0. );
		BOOST_REQUIRE( d1D[1] == 1. );
		BOOST_REQUIRE( d1D[9] == 9. );
	}
	{
		multi::array<double, 1> d1D(begin(multi::index_extension(10)), end(multi::index_extension(10)));
		BOOST_REQUIRE( size(d1D) == 10 );
		BOOST_REQUIRE( d1D[0] == 0. );
		BOOST_REQUIRE( d1D[1] == 1. );
		BOOST_REQUIRE( d1D[9] == 9. );
	}
	{
		multi::array<double, 1> d1D(multi::extensions_t<1>{multi::iextension{10}});
		BOOST_REQUIRE( size(d1D) == 10 );

		d1D.assign(begin(extension(d1D)), end(extension(d1D)));
		BOOST_REQUIRE( d1D[0] == 0. );
		BOOST_REQUIRE( d1D[1] == 1. );
		BOOST_REQUIRE( d1D[9] == 9. );
	}
	 {
		multi::array<double, 1> d1D(multi::extensions_t<1>{multi::iextension{10}});
		d1D.assign(extension(d1D));
		BOOST_REQUIRE( d1D[0] == 0. );
		BOOST_REQUIRE( d1D[1] == 1. );
		BOOST_REQUIRE( d1D[9] == 9. );
	}
}

BOOST_AUTO_TEST_CASE(fill_member) {
	namespace multi = boost::multi;
	multi::array<double, 1> d1D = {1., 2., 3., 4.};
	d1D.fill(42.);

	multi::array<double, 2> d2D = {
		{150., 16., 17., 18., 19.},
		{  5.,  5.,  5.,  5.,  5.},
		{100., 11., 12., 13., 14.},
		{ 50.,  6.,  7.,  8.,  9.}
	};
	BOOST_REQUIRE( d2D.elements().size() == d2D.num_elements() );
	BOOST_REQUIRE( d2D.elements().base() == d2D.base() );
	BOOST_REQUIRE( d2D.elements()[3] == 18. );
	BOOST_REQUIRE( &*d2D.elements().begin() == d2D.data_elements() );
	BOOST_REQUIRE( &*d2D.elements().end() == d2D.data_elements() + d2D.num_elements() );
//	std::fill( d2D.elements().begin(), d2D.elements().end() , 99. );
//	multi::adl_fill_n( d2D.elements().begin(), d2D.elements().size(), 99. );
	d2D.elements().fill(99.);

	BOOST_REQUIRE( d2D[1][1] == 99. );
}

BOOST_AUTO_TEST_CASE(fill) {
	namespace multi = boost::multi;

	multi::array<double, 2> d2D = {
		{150., 16., 17., 18., 19.},
		{  5.,  5.,  5.,  5.,  5.},
		{100., 11., 12., 13., 14.},
		{ 50.,  6.,  7.,  8.,  9.}
	};
	using std::all_of;
	BOOST_REQUIRE( all_of(begin(d2D[1]), end(d2D[1]), [](auto const& e){return e==5.;}) );

	using std::fill;
	fill(d2D[1].begin(), d2D[1].end(), 8.);

	BOOST_REQUIRE( all_of(begin(d2D[1]), end(d2D[1]), [](auto const& e){return e==8.;}) );

	fill(begin(rotated(d2D)[1]), end(rotated(d2D)[1]), 8.);
	BOOST_REQUIRE( all_of(begin(rotated(d2D)[1]), end(rotated(d2D)[1]), [](auto&& e){return e==8.;}) );

	fill(begin((d2D<<1)[1]), end((d2D<<1)[1]), 8.);
	BOOST_REQUIRE( all_of(begin((d2D<<1)[1]), end((d2D<<1)[1]), [](auto&& e){return e==8.;}) );

	std::mt19937 g{r()};
	auto rand = [d=std::normal_distribution<>{}, g = std::mt19937{r()}]() mutable{return d(g);};
	multi::array<double, 2> r2D({5, 5});
	std::for_each(begin(r2D), end(r2D), [&](auto&& e){std::generate(begin(e), end(e), rand);});
}

