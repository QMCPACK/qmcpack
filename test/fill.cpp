// Copyright 2019-2024 Alfredo A. Correa
// Copyright 2024 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/core/lightweight_test.hpp>

#include <boost/multi/array.hpp>  // for array, apply, operator==

#include <algorithm>    // for fill, all_of, transform
#include <cstddef>      // for ptrdiff_t
#include <cstdint>      // for uint64_t
#include <iterator>     // for begin, end, size, next
#include <limits>       // for numeric_limits
#include <numeric>      // for accumulate
#include <random>       // for uniform_int_distribution
#include <type_traits>  // for enable_if_t, is_same_v

namespace {

using fnv1a_size = std::uint64_t;

// from Howard Hinnart hash
auto fnv1a(unsigned char const* first, std::ptrdiff_t len, fnv1a_size hash) noexcept {  // NOLINT(bugprone-easily-swappable-parameters)
	return std::accumulate(
		first, std::next(first, len), hash,
		[prime = 1099511628211U](auto acc, auto elem) { return (acc ^ elem) * prime; }
	);
}
}  // namespace

class fnv1a_t {
	fnv1a_size h_ = 14695981039346656037U;  // offset

 public:
	using result_type = fnv1a_size;
	static constexpr auto(min)() {  // paren for MSVC macros
		return (std::numeric_limits<result_type>::min)();
	}
	static constexpr auto(max)() {
		return (std::numeric_limits<result_type>::max)();
	}

	void operator()(unsigned char const* key, std::ptrdiff_t len) noexcept { h_ = fnv1a(key, len, h_); }

	template<class T, std::enable_if_t<std::is_fundamental_v<T>, int> = 0>  // NOLINT(modernize-use-constraints) TODO(correaa) for C++20
	auto operator()(T const& value) noexcept -> decltype(auto) {
		operator()(&value, sizeof(value));
		return *this;
	}
	//  result_type operator()() && noexcept{return h;}
	auto operator()() const& noexcept { return h_; }
	//  explicit operator result_type() && noexcept {return h;}
	explicit operator result_type() const& noexcept { return h_; }
};

#define BOOST_AUTO_TEST_CASE(CasenamE) /**/

auto main() -> int {  // NOLINT(readability-function-cognitive-complexity,bugprone-exception-escape)
	BOOST_AUTO_TEST_CASE(fill_1d_a) {
		namespace multi = boost::multi;

		multi::array<multi::index, 1> d1D(multi::extensions_t<1>{multi::iextension{10}});
		static_assert(std::is_same_v<std::iterator_traits<decltype(begin(d1D))>::value_type, multi::index>);

		using std::copy;
		copy(begin(extension(d1D)), end(extension(d1D)), begin(d1D));
		BOOST_TEST( d1D[0] == 0 );
		BOOST_TEST( d1D[1] == 1 );
		BOOST_TEST( d1D[9] == 9 );

		d1D.assign(extension(d1D));
		BOOST_TEST( d1D[0] == 0 );
		BOOST_TEST( d1D[1] == 1 );
		BOOST_TEST( d1D[9] == 9 );
	}

	BOOST_AUTO_TEST_CASE(fill_1d_b) {
		namespace multi = boost::multi;

		multi::array<multi::index, 1> d1D(begin(multi::index_extension(10)), end(multi::index_extension(10)));
		BOOST_TEST( size(d1D) == 10 );
		BOOST_TEST( d1D[0] == 0 );
		BOOST_TEST( d1D[1] == 1 );
		BOOST_TEST( d1D[9] == 9 );
	}

	BOOST_AUTO_TEST_CASE(fill_1d_c) {
		namespace multi = boost::multi;

		multi::array<multi::index, 1> d1D(multi::extensions_t<1>{multi::iextension{10}});
		BOOST_TEST( size(d1D) == 10 );

		d1D.assign(begin(extension(d1D)), end(extension(d1D)));
		BOOST_TEST( d1D[0] == 0 );
		BOOST_TEST( d1D[1] == 1 );
		BOOST_TEST( d1D[9] == 9 );
	}

	BOOST_AUTO_TEST_CASE(fill_1d_d) {
		namespace multi = boost::multi;

		multi::array<multi::index, 1> d1D(multi::extensions_t<1>{multi::iextension{10}});
		d1D.assign(extension(d1D));
		BOOST_TEST( d1D[0] == 0 );
		BOOST_TEST( d1D[1] == 1 );
		BOOST_TEST( d1D[9] == 9 );
	}

	BOOST_AUTO_TEST_CASE(fill_member) {
		namespace multi = boost::multi;

		multi::array<int, 1> d1D = {10, 20, 30, 40};
		d1D.fill(420);

		multi::array<int, 2> d2D = {
			{1500, 160, 170, 180, 190},
			{  50,  50,  50,  50,  50},
			{1000, 110, 120, 130, 140},
			{ 500,  60,  70,  80,  90},
		};

		BOOST_TEST(   d2D.elements().size()  == d2D.num_elements()  );
		BOOST_TEST(   d2D.elements().base()  == d2D.base()          );
		BOOST_TEST(   d2D.elements()[3]      == 180                 );

		std::fill(d2D.elements().begin(), d2D.elements().end(), 990);

		BOOST_TEST( d2D[1][1] == 990 );
	}

	BOOST_AUTO_TEST_CASE(simple_fill) {
		namespace multi = boost::multi;

		multi::array<int, 1> d1D = {10, 20, 30, 40};
		std::fill_n(d1D.begin(), d1D.size(), 420);

		multi::array<int, 2> d2D = {
			{1500, 160, 170, 180, 190},
			{  50,  50,  50,  50,  50},
			{1000, 110, 120, 130, 140},
			{ 500,  60,  70,  80,  90},
		};

		BOOST_TEST(   d2D.elements().size()  == d2D.num_elements()  );
		BOOST_TEST(   d2D.elements().base()  == d2D.base()          );
		BOOST_TEST(   d2D.elements()[3]      == 180                 );

		std::fill(d2D.elements().begin(), d2D.elements().end(), 990);

		BOOST_TEST( d2D[1][1] == 990 );
	}

	BOOST_AUTO_TEST_CASE(fill) {
		std::random_device randdev;

		namespace multi = boost::multi;

		multi::array<int, 2> d2D = {
			{1500, 160, 170, 180, 190},
			{  50,  50,  50,  50,  50},
			{1000, 110, 120, 130, 140},
			{ 500,  60,  70,  80,  90},
		};
		using std::all_of;
		BOOST_TEST( all_of(begin(d2D[1]), end(d2D[1]), [](auto const& elem) { return elem == 50;}) );

		using std::fill;
		fill(d2D[1].begin(), d2D[1].end(), 80);

		BOOST_TEST( all_of(begin(d2D[1]), end(d2D[1]), [](auto const& elem) { return elem == 80;}) );

		fill(begin(d2D.rotated()[1]), end(d2D.rotated()[1]), 80);
		BOOST_TEST( all_of(begin(d2D.rotated()[1]), end(d2D.rotated()[1]), [](auto&& elem) { return elem == 80;}) );

		fill(begin((d2D.rotated())[1]), end((d2D.rotated())[1]), 80);
		BOOST_TEST( all_of(begin((d2D.rotated())[1]), end((d2D.rotated())[1]), [](auto&& elem) { return elem == 80;}) );

		auto rand = [gauss = std::uniform_int_distribution<>(0, 10), gen = std::mt19937_64(randdev())]() mutable { return gauss(gen); };  // NOSONAR

		multi::array<int, 2> r2D({5, 5});
		std::for_each(begin(r2D), end(r2D), [&](decltype(r2D)::reference elem) { std::generate(begin(elem), end(elem), rand); });
	}

	BOOST_AUTO_TEST_CASE(fill_1D) {
		namespace multi = boost::multi;

		multi::array<double, 1> const arr = {1.0, 2.0, 3.0};

		multi::array<double, 2> arr2({10, 3});

		std::fill(begin(arr2), end(arr2), arr);

		BOOST_TEST( arr2[0] == arr );
		BOOST_TEST( arr2[1] == arr );

		BOOST_TEST( arr2[9] == arr );
	}

	BOOST_AUTO_TEST_CASE(fill_n_1D) {
		namespace multi = boost::multi;

		multi::array<int, 1> arr = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};

		std::fill_n(arr.dropped(3).begin(), 4, 99);
		// std::fill_n(arr.begin() + 3, 4, 99);

		BOOST_TEST( arr[2] ==  2 );
		BOOST_TEST( arr[3] == 99 );
		BOOST_TEST( arr[4] == 99 );
		BOOST_TEST( arr[5] == 99 );
		BOOST_TEST( arr[6] == 99 );
		BOOST_TEST( arr[7] ==  7 );
	}

	return boost::report_errors();
}
