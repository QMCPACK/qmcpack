// Copyright 2018-2024 Alfredo A. Correa
// Copyright 2024 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/multi/array_ref.hpp>  // for array_ptr, array_ref, subarray

#include <array>        // for array
#include <iterator>     // for iterator_traits
#include <memory>       // for allocator
#include <type_traits>  // for is_same, is_convertible, enable...

namespace multi = boost::multi;

namespace minimalistic {

template<class T>
class ptr : public std::iterator_traits<T*> {  // minimalistic pointer
	using underlying_type = T*;
	underlying_type impl_;
	template<class> friend class ptr;

 public:
	ptr() = default;
	constexpr explicit ptr(T* impl) : impl_{ impl } {}
	template<class U
		, class = std::enable_if_t<std::is_convertible_v<U*, T*>>  // NOLINT(modernize-use-constraints) TODO(correaa) for C++20
	>
	// cppcheck-suppress [noExplicitConstructor,unmatchedSuppression]
	ptr(ptr<U> const& other) : impl_{ other.impl_ } {}  //  NOLINT(google-explicit-constructor, hicpp-explicit-conversions)  // NOSONAR(cpp:S1709)
	using typename std::iterator_traits<T*>::reference;
	using typename std::iterator_traits<T*>::difference_type;

	// NOLINTNEXTLINE(fuchsia-overloaded-operator, fuchsia-trailing-return): operator* used because this class simulates a pointer, trailing return helps
	constexpr auto operator*() const -> reference { return *impl_; }

	#if defined(__clang__)
	#pragma clang diagnostic push
	#pragma clang diagnostic ignored "-Wunknown-warning-option"
	#pragma clang diagnostic ignored "-Wunsafe-buffer-usage"
	#endif

	// NOLINTNEXTLINE(fuchsia-overloaded-operator, cppcoreguidelines-pro-bounds-pointer-arithmetic): operator+ is overloaded to simulate a pointer
	constexpr auto operator+(difference_type n) const { return ptr{ impl_ + n }; }
	// NOLINTNEXTLINE(fuchsia-overloaded-operator, cppcoreguidelines-pro-bounds-pointer-arithmetic): operator+ is overloaded to simulate a pointer
	constexpr auto operator-(difference_type n) const { return ptr{ impl_ - n }; }

	#if defined(__clang__)
	#pragma clang diagnostic pop
	#endif

	//  T& operator[](difference_type n) const{return impl_[n];} // optional
	using default_allocator_type = std::allocator<T>;

	template<class T2> auto operator==(ptr<T2> const& other) const& { return impl_ == other.impl_; }
	template<class> friend class ptr2;
};

template<class T>
class ptr2 : public std::iterator_traits<T*> {  // minimalistic pointer
	T* impl_;

 public:
	constexpr explicit ptr2(T* impl) : impl_{ impl } {}
	constexpr explicit ptr2(ptr<T> const& other) : impl_{ other.impl_ } {}
	template<class U, class = std::enable_if_t<std::is_convertible_v<U*, T*>>>  // NOLINT(modernize-use-constraints) TODO(correaa) for C++20
	// cppcheck-suppress [noExplicitConstructor, unmatchedSuppression]
	ptr2(ptr2<U> const& other) : impl_{ other.impl_ } {}  // NOLINT(google-explicit-constructor, hicpp-explicit-conversions)  // NOSONAR(cpp:S1709)

	using typename std::iterator_traits<T*>::reference;
	using typename std::iterator_traits<T*>::difference_type;

	// NOLINTNEXTLINE(fuchsia-overloaded-operator, fuchsia-trailing-return): operator* used because this class simulates a pointer, trailing return helps
	constexpr auto operator*() const -> reference { return *impl_; }

	#if defined(__clang__)
	#pragma clang diagnostic push
	#pragma clang diagnostic ignored "-Wunknown-warning-option"
	#pragma clang diagnostic ignored "-Wunsafe-buffer-usage"
	#endif

	// NOLINTNEXTLINE(fuchsia-overloaded-operator, cppcoreguidelines-pro-bounds-pointer-arithmetic): operator+ is overloaded to simulate a pointer
	constexpr auto operator+(difference_type n) const { return ptr2{ impl_ + n }; }
	// NOLINTNEXTLINE(fuchsia-overloaded-operator, cppcoreguidelines-pro-bounds-pointer-arithmetic): operator+ is overloaded to simulate a pointer
	constexpr auto operator-(difference_type n) const { return ptr2{ impl_ - n }; }

	#if defined(__clang__)
	#pragma clang diagnostic pop
	#endif

	//  T& operator[](std::ptrdiff_t n) const{return impl_[n];}  // optional
	using default_allocator_type = std::allocator<T>;
};

}  // end namespace minimalistic

#include <boost/core/lightweight_test.hpp>
#define BOOST_AUTO_TEST_CASE(CasenamE) /**/

auto main() -> int {  // NOLINT(readability-function-cognitive-complexity,bugprone-exception-escape)
BOOST_AUTO_TEST_CASE(test_minimalistic_ptr) {
	std::array<double, 400> buffer{};
	BOOST_TEST( buffer.size() == 400 );

	using pointer_type = minimalistic::ptr<double>;
	multi::array_ptr<double, 2, pointer_type> const CCP(pointer_type{ buffer.data() }, { 20, 20 });
	(*CCP)[2];  // requires operator+
	(*CCP)[1][1];
	(*CCP)[1][1] = 9;
	BOOST_TEST( &(*CCP)[1][1] == &buffer[21] );

	// auto&& CC2 = (*CCP).static_array_cast<double, minimalistic::ptr2<double>>();
	auto&& CC2 = CCP->static_array_cast<double, minimalistic::ptr2<double>>();
	BOOST_TEST( &CC2[1][1] == &(*CCP)[1][1] );

	static_assert(std::is_convertible<double*, double const*>{}, "!");

	minimalistic::ptr<double> const       pd{ nullptr };
	minimalistic::ptr<double const> const pcd = pd;
	BOOST_TEST( pcd == pd );

	{
		auto&& REF = *CCP;
		(void)REF;
		static_assert( std::is_same_v<decltype(REF.partitioned(2).partitioned(2).base()), minimalistic::ptr<double>> );
	}
	{
		auto const& REF = *CCP;
		(void)REF;
		static_assert( std::is_same_v<decltype(REF.partitioned(2).partitioned(2).base()), minimalistic::ptr<double const>> );
	}
}
return boost::report_errors();}
