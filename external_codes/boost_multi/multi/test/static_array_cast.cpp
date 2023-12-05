// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2019-2022 Alfredo A. Correa

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi static array cast"
#include<boost/test/unit_test.hpp>

#include "multi/array.hpp"
#include "multi/config/NO_UNIQUE_ADDRESS.hpp"

#include<numeric>

namespace multi = boost::multi;

template<class It, class F> class involuter;

template<class Ref, class Involution>
class involuted {
	Ref r_;
	MULTI_NO_UNIQUE_ADDRESS Involution f_;

 public:
	using decay_type = std::decay_t<decltype(std::declval<Involution>()(std::declval<Ref>()))>;

	constexpr involuted(Ref ref, Involution fun) : r_{std::forward<Ref>(ref)}, f_{fun} {}
	constexpr explicit involuted(Ref ref) : r_{std::forward<Ref>(ref)}, f_{} {}
	involuted(involuted const&) = default;
	involuted(involuted&&) noexcept = default;
	constexpr auto operator=(involuted const& other) = delete;
	~involuted() = default;
	// NOLINTNEXTLINE(google-explicit-constructor,hicpp-explicit-conversions): simulates a reference
	constexpr operator decay_type() const& {return f_(r_);}
	// NOLINTNEXTLINE(google-runtime-operator,fuchsia-overloaded-operator): simulates reference
	constexpr auto operator&() && -> decltype(auto) {return involuter<decltype(&std::declval<Ref>()), Involution>{&r_, f_};} //  NOLINT(runtime/operator)
	// NOLINTNEXTLINE(fuchsia-trailing-return,-warnings-as-errors): trailing return helps reading
	template<class DecayType> constexpr auto operator=(DecayType&& other) & -> involuted& {r_ = f_(std::forward<DecayType>(other)); return *this;}
	// NOLINTNEXTLINE(fuchsia-trailing-return): trailing return helps reading
	constexpr auto operator=(involuted&& other)& noexcept -> involuted& = default;

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
	It it_;
	MULTI_NO_UNIQUE_ADDRESS F f_;
	template<class, class> friend class involuter;
//	template<class From, std::enable_if_t<std::is_convertible<From, It>{}, int> =0>
//	static constexpr auto implicit_cast(From&& f) {return static_cast<It>(f);}

 public:
	using pointer           = involuter<typename std::iterator_traits<It>::pointer, F>;
	using element_type      = typename std::pointer_traits<It>::element_type;
	using difference_type   = typename std::pointer_traits<It>::difference_type;
	template<class U>
	using rebind            = involuter<typename std::pointer_traits<It>::template rebind<U>, F>;
	using reference         = involuted<typename std::iterator_traits<It>::reference, F>;
	using value_type        = typename std::iterator_traits<It>::value_type;
	using iterator_category = typename std::iterator_traits<It>::iterator_category;
	explicit constexpr involuter(It it) : it_{std::move(it)}, f_{} {}  // NOLINT(readability-identifier-length) clang-tidy 14 bug
	constexpr involuter(It it, F fun) : it_{std::move(it)}, f_{std::move(fun)} {}
//	involuter(involuter const& other) = default;
	// NOLINTNEXTLINE(google-explicit-constructor, hicpp-explicit-conversions): this is needed to make involuter<T> implicitly convertible to involuter<T const>
	template<class Other> constexpr involuter(involuter<Other, F> const& other) : it_{multi::implicit_cast<It>(other.it_)}, f_{other.f_} {}
//	auto operator=(involuter const& other) -> involuter& = default;
	constexpr auto operator*() const {return reference{*it_, f_};}
	constexpr auto operator==(involuter const& other) const {return it_ == other.it_;}
	constexpr auto operator!=(involuter const& other) const {return it_ != other.it_;}
	constexpr auto operator+=(typename involuter::difference_type n) -> decltype(auto) {it_+=n; return *this;}
	constexpr auto operator+ (typename involuter::difference_type n) const {return involuter{it_+n, f_};}
	constexpr auto operator- (typename involuter::difference_type n) const {return involuter{it_-n, f_};}
	constexpr auto operator-(involuter const& other) const {return it_ - other.it_;}
	constexpr auto operator->() const {return pointer{&*it_, f_};}
//	~involuter() = default;
	constexpr auto operator[](typename involuter::difference_type n) const {return reference{*(it_ + n), f_};}
};

#if defined(__cpp_deduction_guides)
template<class T, class F> involuted(T&&, F)->involuted<T const, F>;
#endif

template<class Ref> using negated = involuted<Ref, std::negate<>>;
template<class It>  using negater = involuter<It, std::negate<>>;

BOOST_AUTO_TEST_CASE(multi_array_involution) {
	double doub = 5;

	auto&& cee = involuted<double&, std::negate<>>{doub};
	BOOST_REQUIRE( cee == -5. );

	cee = 10.;
	BOOST_REQUIRE( doub = -10. );

	auto m5 = involuted<double, std::negate<>>(5.);
	BOOST_REQUIRE( m5 == -5. );
}

BOOST_AUTO_TEST_CASE(static_array_cast) {
	multi::static_array<double, 1> arr = { 0.,  1.,  2.,  3.,  4.};
	auto&& ref = arr.static_array_cast<double, double const*>();
	BOOST_REQUIRE( &ref[2] == &arr    [2] );
	BOOST_REQUIRE( &arr    [2] == &ref[2] );

	BOOST_REQUIRE( std::equal(begin(ref), end(ref), begin(arr), end(arr)) );
	BOOST_REQUIRE( ref == arr     );
	BOOST_REQUIRE( arr == ref );
}

BOOST_AUTO_TEST_CASE(static_array_cast_2) {
	multi::array<double, 2> arr({2, 5});
	std::iota(arr.elements().begin(), arr.elements().end(), 0.);

	auto&& ref = arr.static_array_cast<double, double const*>();
	BOOST_REQUIRE( ref[1][1] == arr[1][1] );
	BOOST_REQUIRE( std::equal(begin(ref[1]), end(ref[1]), begin(arr[1]), end(arr[1])) );
	BOOST_REQUIRE( ref[1] == arr[1] );
	BOOST_REQUIRE( std::equal(begin(ref), end(ref), begin(arr), end(arr)) );
	BOOST_REQUIRE( ref == arr     );
	BOOST_REQUIRE( arr     == ref );
}

BOOST_AUTO_TEST_CASE(static_array_cast_3) {
{
	multi::static_array<double, 1> arr  { {  0.,   1.,   2.,   3.,  4.} };
	multi::static_array<double, 1> arr2 = { -0.,  -1.,  -2.,  -3., -4.};
	auto&& neg_arr = multi::static_array_cast<double, involuter<double*, std::negate<>>>(arr);
	BOOST_REQUIRE( neg_arr[2] == arr2[2] );
	BOOST_REQUIRE( arr2[2] == neg_arr[2] );
	BOOST_REQUIRE( std::equal(begin(neg_arr), end(neg_arr), begin(arr2), end(arr2)) );
	BOOST_REQUIRE( neg_arr == arr2 );
	BOOST_REQUIRE( arr2 == neg_arr );
}
{
	multi::static_array<double, 2> arr({4, 5}, 0.);
	std::iota(elements(arr).begin(), elements(arr).end(), 0.);

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
