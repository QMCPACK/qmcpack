// Copyright 2019-2024 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/multi/array.hpp>
#include <boost/multi/detail/config/NO_UNIQUE_ADDRESS.hpp>  // TODO(correaa) remove in c++20

#include <algorithm>    // for equal
#include <cassert>      // for assert
#include <functional>   // for negate  // IWYU pragma: keep
#include <iterator>     // for begin, end
#include <memory>       // for pointer_t...
#include <numeric>      // for iota
#include <type_traits>  // for decay_t
#include <utility>      // for move, dec...

namespace multi = boost::multi;

template<class Ref, class Involution>
class involuted {
	Ref                                      r_;  // NOLINT(cppcoreguidelines-avoid-const-or-ref-data-members)
	BOOST_MULTI_NO_UNIQUE_ADDRESS Involution f_;  // TODO(correaa) put nounique members first?

 public:
	using decay_type = std::decay_t<decltype(std::declval<Involution>()(std::declval<Ref>()))>;

	constexpr involuted(Ref ref, Involution fun) : r_{ref}, f_{fun} {}
	constexpr explicit involuted(Ref ref) : r_{ref}, f_{} {}

	involuted(involuted const&)     = default;
	involuted(involuted&&) noexcept = default;

	constexpr auto operator=(involuted const& other) = delete;

	~involuted() = default;

	// NOLINTNEXTLINE(google-explicit-constructor,hicpp-explicit-conversions)
	constexpr operator decay_type() const& noexcept { return f_(r_); }  // NOSONAR(cpp:S1709) simulates a reference

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
	It                              it_;
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
	template<class Other> constexpr involuter(involuter<Other, F> const& other)  // NOSONAR(cpp:S1709)
	: it_{multi::detail::implicit_cast<It>(other.it_)}, f_{other.f_} {}

	constexpr auto operator*() const { return reference{*it_, f_}; }
	constexpr auto operator->() const { return pointer{&*it_, f_}; }

	constexpr auto operator==(involuter const& other) const { return it_ == other.it_; }
	constexpr auto operator!=(involuter const& other) const { return it_ != other.it_; }
	constexpr auto operator<(involuter const& other) const { return it_ < other.it_; }

	#if defined(__clang__)
	#pragma clang diagnostic push
	#pragma clang diagnostic ignored "-Wunknown-warning-option"
	#pragma clang diagnostic ignored "-Wunsafe-buffer-usage"
	#endif

	constexpr auto operator+=(typename involuter::difference_type n) -> decltype(auto) {
		it_ += n;
		return *this;
	}

	#if defined(__clang__)
	#pragma clang diagnostic pop
	#endif

	#if defined(__clang__)
	#pragma clang diagnostic push
	#pragma clang diagnostic ignored "-Wunknown-warning-option"
	#pragma clang diagnostic ignored "-Wunsafe-buffer-usage"
	#endif

	constexpr auto operator+(typename involuter::difference_type n) const { return involuter{it_ + n, f_}; }
	constexpr auto operator-(typename involuter::difference_type n) const { return involuter{it_ - n, f_}; }

	constexpr auto operator[](typename involuter::difference_type n) const { return reference{*(it_ + n), f_}; }

	#if defined(__clang__)
	#pragma clang diagnostic pop
	#endif

	constexpr auto operator-(involuter const& other) const { return it_ - other.it_; }
};

#if defined(__cpp_deduction_guides)
template<class T, class F> involuted(T&&, F) -> involuted<T const, F>;
#endif

template<class Ref> using negated = involuted<Ref, std::negate<>>;
template<class Ptr> using negater = involuter<Ptr, std::negate<>>;

#include <boost/core/lightweight_test.hpp>
#define BOOST_AUTO_TEST_CASE(CasenamE) /**/

auto main() -> int {  // NOLINT(readability-function-cognitive-complexity,bugprone-exception-escape)
BOOST_AUTO_TEST_CASE(multi_array_involution) {
	int doub = 50;

	auto&& cee = involuted<int&, std::negate<>>{doub};
	BOOST_TEST( cee == -50 );

	cee = 100;
	BOOST_TEST( doub == -100 );

	auto m5 = involuted<int, std::negate<>>(50);
	BOOST_TEST( m5 == -50 );
}

BOOST_AUTO_TEST_CASE(static_array_cast) {
	multi::static_array<double, 1> arr = {0.0, 1.0, 2.0, 3.0, 4.0};

	auto&& ref = arr.static_array_cast<double, double const*>();

	BOOST_TEST( &ref[2] == &arr[2] );
	BOOST_TEST( &arr[2] == &ref[2] );

	BOOST_TEST( std::equal(begin(ref), end(ref), begin(arr), end(arr)) );

	BOOST_TEST( ref == arr() );
	BOOST_TEST( arr() == ref );

	BOOST_TEST( ref == arr );
	BOOST_TEST( arr == ref );
}

BOOST_AUTO_TEST_CASE(static_array_cast_2) {
	multi::array<int, 2> arr({2, 5});
	std::iota(arr.elements().begin(), arr.elements().end(), 0);

	auto&& ref = arr.static_array_cast<int, int const*>();

	BOOST_TEST( ref[1][1] == arr[1][1] );
	BOOST_TEST( std::equal(begin(ref[1]), end(ref[1]), begin(arr[1]), end(arr[1])) );
	BOOST_TEST( ref[1] == arr[1] );

	BOOST_TEST( std::equal(begin(ref), end(ref), begin(arr), end(arr)) );

	BOOST_TEST( ref == arr );
	BOOST_TEST( arr == ref );
}

BOOST_AUTO_TEST_CASE(static_array_cast_3) {
	{
		multi::static_array<int, 1> const arr  = {+00, +10, +20, +30, +40};
		multi::static_array<int, 1>       arr2 = {-00, -10, -20, -30, -40};

		auto&& neg_arr = multi::static_array_cast<int, involuter<int*, std::negate<>>>(arr);

		BOOST_TEST( neg_arr[2] == arr2[2] );
		BOOST_TEST( arr2[2] == neg_arr[2] );
		BOOST_TEST( std::equal(begin(neg_arr), end(neg_arr), begin(arr2), end(arr2)) );
		BOOST_TEST( neg_arr == arr2 );
		BOOST_TEST( arr2 == neg_arr );
	}
	{
		multi::static_array<int, 2> arr({4, 5}, 0);
		std::iota(elements(arr).begin(), elements(arr).end(), 0);

		multi::array<int, 2> arr2({4, 5});
		std::transform(begin(elements(arr)), end(elements(arr)), begin(elements(arr2)), std::negate<>{});

		auto&& neg_arr = arr.static_array_cast<int, negater<int*>>();

		BOOST_TEST( neg_arr[1][1] == arr2[1][1] );
		BOOST_TEST( arr2[1][1] == neg_arr[1][1] );

		BOOST_TEST( std::equal(begin(arr2[1]), end(arr2[1]), begin(neg_arr[1]), end(neg_arr[1])) );

		BOOST_TEST( arr2[1] == neg_arr[1] );
		BOOST_TEST( neg_arr[1] == arr2[1] );

		BOOST_TEST( std::equal(begin(arr2), end(arr2), begin(neg_arr), end(neg_arr)) );
		BOOST_TEST( neg_arr == arr2 );
		BOOST_TEST( arr2 == neg_arr );
	}
}
return boost::report_errors();}
