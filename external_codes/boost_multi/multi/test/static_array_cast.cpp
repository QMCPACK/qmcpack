// Copyright 2019-2025 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/multi/array.hpp>
#include <boost/multi/detail/config/NO_UNIQUE_ADDRESS.hpp>  // TODO(correaa) remove in c++20

#include <boost/core/lightweight_test.hpp>

#include <algorithm>    // for equal
#include <array>        // for array
#include <cassert>      // for assert
#include <functional>   // for negate  // IWYU pragma: keep
#include <iterator>     // for begin, end
#include <memory>       // for pointer_t...
#include <numeric>      // for iota
#include <type_traits>  // for decay_t
#include <utility>      // for move, dec...

namespace multi = boost::multi;

template<class Ref, class Involution>
class involuted {  // NOLINT(misc-use-internal-linkage)
#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable : 4820)  // '3' bytes padding added after data member 'involuted<int,std::negate<void>>::f_'
#endif
	BOOST_MULTI_NO_UNIQUE_ADDRESS Involution f_;  // TODO(correaa) put nounique members first?
#ifdef _MSC_VER
#pragma warning(pop)
#endif
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wpadded"
#endif
	Ref r_;  // NOLINT(cppcoreguidelines-avoid-const-or-ref-data-members)
#ifdef __clang__
#pragma clang diagnostic pop
#endif

 public:
	using decay_type = std::decay_t<decltype(std::declval<Involution>()(std::declval<Ref>()))>;

	constexpr involuted(Ref ref, Involution fun)
	// : r_{ref}, f_{fun} {}
	: f_{fun}, r_{ref} {}

	constexpr explicit involuted(Ref ref)
	// : r_{ref}, f_{} {}
	: f_{}, r_{ref} {}

	involuted(involuted const&)     = default;
	involuted(involuted&&) noexcept = default;

	constexpr auto operator=(involuted const&) = delete;
	constexpr auto operator=(involuted&&)      = delete;

	~involuted() = default;

	// NOLINTNEXTLINE(google-explicit-constructor,hicpp-explicit-conversions)
	constexpr operator decay_type() const& noexcept { return f_(r_); }  // NOSONAR(cpp:S1709) simulates a reference

	// NOLINTNEXTLINE(fuchsia-trailing-return,-warnings-as-errors): trailing return helps reading
	template<class DecayType> constexpr auto operator=(DecayType&& other) & -> involuted& {
		r_ = f_(std::forward<DecayType>(other));
		return *this;
	}

	friend auto operator==(involuted const& self, involuted const& other) -> bool {
		assert(self.f_ == other.f_);
		return self.r_ == other.r_;
	}
	friend auto operator!=(involuted const& self, involuted const& other) -> bool {
		assert(self.f_ == other.f_);
		return self.r_ != other.r_;
	}
};

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wpadded"
#endif
template<class It, class F>
class involuter {  // NOLINT(misc-use-internal-linkage)
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

	// vvv this is needed to make involuter<T> implicitly convertible to involuter<T const>
	// cppcheck-suppress noExplicitConstructor ;  // NOLINTNEXTLINE(google-explicit-constructor, hicpp-explicit-conversions)
	template<class Other> constexpr involuter(involuter<Other, F> const& other)  // NOSONAR(cpp:S1709)
	: it_{multi::detail::implicit_cast<It>(other.it_)}, f_{other.f_} {}

	constexpr auto operator*() const { return reference{*it_, f_}; }
	constexpr auto operator->() const { return pointer{&*it_, f_}; }  // cppcheck-suppress redundantPointerOp ; lib idiom

	constexpr auto operator==(involuter const& other) const { return it_ == other.it_; }
	constexpr auto operator!=(involuter const& other) const { return it_ != other.it_; }
	constexpr auto operator<(involuter const& other) const { return it_ < other.it_; }

#if defined(__clang__) && (__clang_major__ >= 16) && !defined(__INTEL_LLVM_COMPILER)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunsafe-buffer-usage"
#endif

	constexpr auto operator+=(typename involuter::difference_type n) -> decltype(auto) {
		it_ += n;  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
		return *this;
	}

	// NOLINTBEGIN(cppcoreguidelines-pro-bounds-pointer-arithmetic)
	constexpr auto operator+(typename involuter::difference_type n) const { return involuter{it_ + n, f_}; }
	constexpr auto operator-(typename involuter::difference_type n) const { return involuter{it_ - n, f_}; }

	constexpr auto operator[](typename involuter::difference_type n) const { return reference{*(it_ + n), f_}; }
	// NOLINTEND(cppcoreguidelines-pro-bounds-pointer-arithmetic)

	friend constexpr auto operator+(typename involuter::difference_type n, involuter const& self) { return self + n; }

#if defined(__clang__) && (__clang_major__ >= 16) && !defined(__INTEL_LLVM_COMPILER)
#pragma clang diagnostic pop
#endif

	constexpr auto operator-(involuter const& other) const { return it_ - other.it_; }
};
#ifdef __clang__
#pragma clang diagnostic pop
#endif

#ifdef __cpp_deduction_guides
template<class T, class F> involuted(T&&, F) -> involuted<T const, F>;  // NOLINT(misc-use-internal-linkage) bug in clang-tidy 19
#endif

template<class Ref> using negated = involuted<Ref, std::negate<>>;
template<class Ptr> using negater = involuter<Ptr, std::negate<>>;

auto main() -> int {  // NOLINT(readability-function-cognitive-complexity,bugprone-exception-escape)
	// BOOST_AUTO_TEST_CASE(multi_array_involution)
	{
		int doub = 50;

		auto&& cee = involuted<int&, std::negate<>>{doub};
		BOOST_TEST( cee == -50 );

		cee = 100;
		BOOST_TEST( doub == -100 );

		auto m5 = involuted<int, std::negate<>>(50);
		BOOST_TEST( m5 == -50 );
	}

	// BOOST_AUTO_TEST_CASE(dynamic_array_cast)
	{
		multi::dynamic_array<double, 1> arr = {0.0, 1.0, 2.0, 3.0, 4.0};

		auto&& ref = arr.static_array_cast<double, double const*>();

		BOOST_TEST( &ref[2] == &arr[2] );
		BOOST_TEST( &arr[2] == &ref[2] );

		BOOST_TEST( std::equal(begin(ref), end(ref), begin(arr), end(arr)) );

		BOOST_TEST( ref == arr() );
		BOOST_TEST( arr() == ref );

		BOOST_TEST( ref == arr );
		BOOST_TEST( arr == ref );
	}

	// BOOST_AUTO_TEST_CASE(static_array_cast_2)
	{
		multi::array<int, 2> arr({2, 5});
		std::iota(arr.elements().begin(), arr.elements().end(), 0);

		auto&& ref = arr.static_array_cast<int, int const*>();

		BOOST_TEST( ref[1][1] == arr[1][1] );
		BOOST_TEST( std::equal(ref[1].begin(), ref[1].end(), arr[1].begin(), arr[1].end()) );
		BOOST_TEST( ref[1] == arr[1] );

#if !defined(_MSC_VER) && !defined(__NVCC__)
		BOOST_TEST( std::equal(ref.begin(), ref.end(), arr.begin(), arr.end()) );  // NOLINT(llvm-use-ranges,modernize-use-ranges) for C++20
#endif
		// ^^^ this doesn't work on MSVC+NVCC in C++20 because it tries to generate this type:
		// using coty = std::common_reference<
		// 	boost::multi::subarray<int, 1LL, const int *, boost::multi::layout_t<1LL, boost::multi::size_type>> &&,
		//  	boost::multi::array<int, 1LL, std::allocator<int>> &
		// >::type;
		// instantiation of type "std::_Cond_res<boost::multi::subarray<int, 1LL, const int *, boost::multi::layout_t
		// <1LL, boost::multi::size_type>> &&, boost::multi::array<int, 1LL, std::allocator<int>> &>" at line 1405
		// instantiation of class "std::_Common_reference2C<_Ty1, _Ty2> [with _Ty1=boost::multi::subarray
		// <int, 1LL, const int *, boost::multi::layout_t<1LL, boost::multi::size_type>> &&, _Ty2=
		// boost::multi::array<int, 1LL, std::allocator<int>> &]" at line 1414
		// instantiation of class "std::_Common_reference2B<_Ty1, _Ty2> [with _Ty1=boost::multi::subarray
		// <int, 1LL, const int *, boost::multi::layout_t<1LL, boost::multi::size_type>> &&, _Ty2=boost::multi::array<int, 1LL, std::allocator<int>> &]" at line 1426
		// instantiation of class "std::_Common_reference2A<_Ty1, _Ty2> [with _Ty1=boost::multi::subarray
		// <int, 1LL, const int *, boost::multi::layout_t<1LL, boost::multi::size_type>> &&, _Ty2=boost::multi::array<int, 1LL, std::allocator<int>> &]" at line 1474
		// instantiation of class "std::common_reference<_Ty1, _Ty2> [with
		// 	_Ty1=boost::multi::subarray<int, 1LL, const int *, boost::multi::layout_t<1LL, boost::multi::size_type>> &&,
		// 	_Ty2=boost::multi::array<int, 1LL, std::allocator<int>> &]" at line 1313 of
		// C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Tools\MSVC\14.44.35207\include\xutility
		// instantiation of "const __nv_bool std::_Is_ranges_random_iter_v [with _Iter=
		// 	boost::multi::array_iterator<int, 2LL, const int *, false, false, ptrdiff_t>]" at line 5563 of
		// C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Tools\MSVC\14.44.35207\include\xutility
		// instantiation of "__nv_bool std::equal(_InIt1, _InIt1, _InIt2, _InIt2, _Pr) [with
		// 		_InIt1=boost::multi::array_iterator<int, 2LL, const int *, false, false, ptrdiff_t>,
		// 		_InIt2=boost::multi::array_iterator<int, 2LL, int *      , false, false, ptrdiff_t>, _Pr=std::equal_to<void>]"
		// at line 5599 of C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Tools\MSVC\14.44.35207\include\xutility
		// instantiation of "__nv_bool std::equal(_InIt1, _InIt1, _InIt2, _InIt2) [with
		// 	_InIt1=boost::multi::array_iterator<int, 2LL, const int *, false, false, ptrdiff_t>,
		// 	_InIt2=boost::multi::array_iterator<int, 2LL, int *.     , false, false, ptrdiff_t>]"
		// at line 193 of C:\Gitlab-Runner\builds\t3_1sV2uA\0\correaa\boost-multi\test\static_array_cast.cpp

		BOOST_TEST( ref == arr );
		BOOST_TEST( arr == ref );
	}

	// // BOOST_AUTO_TEST_CASE(static_array_cast_3)
	{
		multi::dynamic_array<int, 1> const arr  = {+00, +10, +20, +30, +40};
		multi::dynamic_array<int, 1>       arr2 = {-00, -10, -20, -30, -40};

		auto&& neg_arr = multi::static_array_cast<int, involuter<int*, std::negate<>>>(arr);

		BOOST_TEST( neg_arr[2] == arr2[2] );
		BOOST_TEST( arr2[2] == neg_arr[2] );
		BOOST_TEST( std::equal( neg_arr.begin(), neg_arr.end(), arr2.begin(), arr2.end() ) );  // NOLINT(llvm-use-ranges,modernize-use-ranges) for C++20
		BOOST_TEST( neg_arr == arr2 );
		BOOST_TEST( arr2 == neg_arr );
	}
	{
		multi::dynamic_array<int, 2> arr({4, 5}, 0);
		std::iota(elements(arr).begin(), elements(arr).end(), 0);

		multi::array<int, 2> arr2({4, 5});
		std::transform(begin(elements(arr)), end(elements(arr)), begin(elements(arr2)), std::negate<>{});

		auto&& neg_arr = arr.static_array_cast<int, negater<int*>>();

		BOOST_TEST( neg_arr[1][1] == arr2[1][1] );
		BOOST_TEST( arr2[1][1] == neg_arr[1][1] );

		BOOST_TEST( std::equal(begin(arr2[1]), end(arr2[1]), begin(neg_arr[1]), end(neg_arr[1])) );

		BOOST_TEST( arr2[1] == neg_arr[1] );
		BOOST_TEST( neg_arr[1] == arr2[1] );

		BOOST_TEST( std::equal(arr2.begin(), arr2.end(), neg_arr.begin(), neg_arr.end()) );  // NOLINT(llvm-use-ranges,modernize-use-ranges) for C++20
		BOOST_TEST( neg_arr == arr2 );
		BOOST_TEST( arr2 == neg_arr );
	}
	{
		multi::array<int, 2> const arr({3, 4}, multi::uninitialized_elements_t{});
		BOOST_TEST( arr.size() == 3 );
	}
	{
		multi::array<int, 2> const arr({3, 4});
		// std::cout << arr[0][0] << std::endl;  // ok, gives an error in Valgrind "Uninitialized Memory Read"
		// valgrind output:
		// Memory checking results:
		// Uninitialized Memory Conditional - 6
		// Uninitialized Memory Read - 2
		BOOST_TEST( arr.size() == 3 );
	}
	{
		multi::array<int, 2> const arr1({
			{3, 5}
		});  // this array has 1-row
		multi::array<int, 2> const arr2 = {
			{3, 5}
		};  // this array has 1-row

		BOOST_TEST( arr1 == arr2  );
	}
	{
		multi::array<int, 2> const arr(std::array<multi::ssize_t, 2>{
			{3, 4}
		});
		BOOST_TEST( arr.size() == 3 );
	}
	{
		multi::array<int, 2> const arr({3, 4}, multi::uninitialized_elements);
		// std::cout << arr[0][0] << std::endl;  // ok, gives an error in Valgrind "Uninitialized Memory Read"
		// valgrind output:
		// Memory checking results:
		// Uninitialized Memory Conditional - 6
		// Uninitialized Memory Read - 2
		BOOST_TEST( arr.size() == 3 );
	}
	{
		// multi::array<std::string, 2> arr( {3, 4}, multi::uninitialized );  // ok, fails compilation because std::string cannot be uninitialized
		// BOOST_TEST( arr.size() == 3 );
	}

	return boost::report_errors();
}
