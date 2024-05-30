// Copyright 2019-2024 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/multi/array.hpp>
#include <boost/multi/detail/config/NO_UNIQUE_ADDRESS.hpp>  // TODO(correaa) remove in c++20

#include <numeric>

// Suppress warnings from boost.test
#if defined(__clang__)
#  pragma clang diagnostic push
#  pragma clang diagnostic ignored "-Wold-style-cast"
#  pragma clang diagnostic ignored "-Wundef"
#  pragma clang diagnostic ignored "-Wconversion"
#  pragma clang diagnostic ignored "-Wsign-conversion"
#  pragma clang diagnostic ignored "-Wfloat-equal"
#elif defined(__GNUC__)
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wold-style-cast"
#  pragma GCC diagnostic ignored "-Wundef"
#  pragma GCC diagnostic ignored "-Wconversion"
#  pragma GCC diagnostic ignored "-Wsign-conversion"
#  pragma GCC diagnostic ignored "-Wfloat-equal"
#endif

#ifndef BOOST_TEST_MODULE
#  define BOOST_TEST_MAIN
#endif

#include <boost/test/unit_test.hpp>

namespace multi = boost::multi;

template<class It, class F> class involuter;

template<class Ref, class Involution>
class involuted {
	Ref                                r_;  // NOLINT(cppcoreguidelines-avoid-const-or-ref-data-members)
	BOOST_MULTI_NO_UNIQUE_ADDRESS Involution f_;  // TODO(correaa) put nounique members first?

 public:
	using decay_type = std::decay_t<decltype(std::declval<Involution>()(std::declval<Ref>()))>;

	constexpr involuted(Ref ref, Involution fun) : r_{ref}, f_{fun} {}
	constexpr explicit involuted(Ref ref) : r_{ref}, f_{} {}

	involuted(involuted const&)     = default;
	involuted(involuted&&) noexcept = default;

	constexpr auto operator=(involuted const& other) = delete;

	~involuted() = default;

	constexpr operator decay_type() const& noexcept { return f_(r_); }  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions)  // NOSONAR(cpp:S1709) simulates a reference
	// NOLINTNEXTLINE(fuchsia-trailing-return,-warnings-as-errors): trailing return helps reading
	template<class DecayType> constexpr auto operator=(DecayType&& other) & -> involuted& {
		r_ = f_(std::forward<DecayType>(other));
		return *this;
	}
	// NOLINTNEXTLINE(fuchsia-trailing-return): trailing return helps reading
	constexpr auto operator=(involuted&& other) & noexcept -> involuted& = default;

	friend auto operator==(involuted const& self, involuted const& other) -> bool {
		assert(self.f_ == other.f_);
		return self.r_ == other.r_;
	}
	friend auto operator!=(involuted const& self, involuted const& other) -> bool {
		assert(self.f_ == other.f_);
		return self.r_ != other.r_;
	}
};

template<class It, class F>
class involuter {
	It                        it_;
	BOOST_MULTI_NO_UNIQUE_ADDRESS F f_;
	template<class, class> friend class involuter;

 public:
	using pointer         = involuter<typename std::iterator_traits<It>::pointer, F>;
	using element_type    = typename std::pointer_traits<It>::element_type;
	using difference_type = typename std::pointer_traits<It>::difference_type;
	template<class U>
	using rebind            = involuter<typename std::pointer_traits<It>::template rebind<U>, F>;
	using reference         = involuted<typename std::iterator_traits<It>::reference, F>;
	using value_type        = typename std::iterator_traits<It>::value_type;
	using iterator_category = typename std::iterator_traits<It>::iterator_category;
	explicit constexpr involuter(It it) : it_{std::move(it)}, f_{} {}  // NOLINT(readability-identifier-length) clang-tidy 14 bug
	constexpr involuter(It it, F fun) : it_{std::move(it)}, f_{std::move(fun)} {}

	// NOLINTNEXTLINE(google-explicit-constructor, hicpp-explicit-conversions): this is needed to make involuter<T> implicitly convertible to involuter<T const>
	template<class Other> constexpr involuter(involuter<Other, F> const& other) : it_{multi::detail::implicit_cast<It>(other.it_)}, f_{other.f_} {}  // NOSONAR(cpp:S1709)

	constexpr auto operator*() const { return reference{*it_, f_}; }
	constexpr auto operator->() const { return pointer{&*it_, f_}; }

	constexpr auto operator==(involuter const& other) const { return it_ == other.it_; }
	constexpr auto operator!=(involuter const& other) const { return it_ != other.it_; }
	constexpr auto operator<(involuter const& other) const { return it_ < other.it_; }

	constexpr auto operator+=(typename involuter::difference_type n) -> decltype(auto) {
		it_ += n;
		return *this;
	}
	constexpr auto operator+(typename involuter::difference_type n) const { return involuter{it_ + n, f_}; }
	constexpr auto operator-(typename involuter::difference_type n) const { return involuter{it_ - n, f_}; }
	constexpr auto operator-(involuter const& other) const { return it_ - other.it_; }

	constexpr auto operator[](typename involuter::difference_type n) const { return reference{*(it_ + n), f_}; }
};

#if defined(__cpp_deduction_guides)
template<class T, class F> involuted(T&&, F) -> involuted<T const, F>;
#endif

template<class Ref> using negated = involuted<Ref, std::negate<>>;
template<class Ptr> using negater = involuter<Ptr, std::negate<>>;

BOOST_AUTO_TEST_CASE(multi_array_involution) {
	double doub = 5;

	auto&& cee = involuted<double&, std::negate<>>{doub};
	BOOST_REQUIRE( cee == -5.0 );

	cee = 10.;
	BOOST_REQUIRE( doub = -10.0 );

	auto m5 = involuted<double, std::negate<>>(5.0);
	BOOST_REQUIRE( m5 == -5.0 );
}

BOOST_AUTO_TEST_CASE(static_array_cast) {
	multi::static_array<double, 1> arr = {0.0, 1.0, 2.0, 3.0, 4.0};

	auto&& ref = arr.static_array_cast<double, double const*>();

	BOOST_REQUIRE( &ref[2] == &arr[2] );
	BOOST_REQUIRE( &arr[2] == &ref[2] );

	BOOST_REQUIRE( std::equal(begin(ref), end(ref), begin(arr), end(arr)) );

	BOOST_REQUIRE( ref == arr() );
	BOOST_REQUIRE( arr() == ref );

	BOOST_REQUIRE( ref == arr );
	BOOST_REQUIRE( arr == ref );
}

BOOST_AUTO_TEST_CASE(static_array_cast_2) {
	multi::array<double, 2> arr({2, 5});
	std::iota(arr.elements().begin(), arr.elements().end(), 0.0);

	auto&& ref = arr.static_array_cast<double, double const*>();

	BOOST_REQUIRE( ref[1][1] == arr[1][1] );
	BOOST_REQUIRE( std::equal(begin(ref[1]), end(ref[1]), begin(arr[1]), end(arr[1])) );
	BOOST_REQUIRE( ref[1] == arr[1] );

	BOOST_REQUIRE( std::equal(begin(ref), end(ref), begin(arr), end(arr)) );

	BOOST_REQUIRE( ref == arr );
	BOOST_REQUIRE( arr == ref );
}

BOOST_AUTO_TEST_CASE(static_array_cast_3) {
	{
		multi::static_array<double, 1> const arr  = {+0.0, +1.0, +2.0, +3.0, +4.0};
		multi::static_array<double, 1>       arr2 = {-0.0, -1.0, -2.0, -3.0, -4.0};

		auto&& neg_arr = multi::static_array_cast<double, involuter<double*, std::negate<>>>(arr);

		BOOST_REQUIRE( neg_arr[2] == arr2[2] );
		BOOST_REQUIRE( arr2[2] == neg_arr[2] );
		BOOST_REQUIRE( std::equal(begin(neg_arr), end(neg_arr), begin(arr2), end(arr2)) );
		BOOST_REQUIRE( neg_arr == arr2 );
		BOOST_REQUIRE( arr2 == neg_arr );
	}
	{
		multi::static_array<double, 2> arr({4, 5}, 0.0);
		std::iota(elements(arr).begin(), elements(arr).end(), 0.0);

		multi::array<double, 2> arr2({4, 5});
		std::transform(begin(elements(arr)), end(elements(arr)), begin(elements(arr2)), std::negate<>{});

		auto&& neg_arr = arr.static_array_cast<double, negater<double*>>();

		BOOST_REQUIRE( neg_arr[1][1] == arr2[1][1] );
		BOOST_REQUIRE( arr2[1][1] == neg_arr[1][1] );

		BOOST_REQUIRE( std::equal(begin(arr2[1]), end(arr2[1]), begin(neg_arr[1]), end(neg_arr[1])) );

		BOOST_REQUIRE( arr2[1] == neg_arr[1] );
		BOOST_REQUIRE( neg_arr[1] == arr2[1] );

		BOOST_REQUIRE( std::equal(begin(arr2), end(arr2), begin(neg_arr), end(neg_arr)) );
		BOOST_REQUIRE( neg_arr == arr2 );
		BOOST_REQUIRE( arr2 == neg_arr );
	}
}
