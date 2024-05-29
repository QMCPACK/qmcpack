// Copyright 2018-2024 Alfredo A. Correa
// Copyright 2024 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/multi/array.hpp>

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
#elif defined(_MSC_VER)
#  pragma warning(push)
#  pragma warning(disable : 4285)  // Recursive return type for fancy_ptr if infix notationis applied
#endif

#ifndef BOOST_TEST_MODULE
#  define BOOST_TEST_MAIN
#endif

#include <boost/test/unit_test.hpp>

#if defined(__clang__)
#  pragma clang diagnostic pop
#elif defined(__GNUC__)
#  pragma GCC diagnostic pop
#endif

namespace fancy {

template<class T> class ref;

template<class T = void> class ptr {  // NOLINT(cppcoreguidelines-special-member-functions,hicpp-special-member-functions)
	static double const value;

 public:
	using difference_type   = std::ptrdiff_t;
	using value_type        = std::decay_t<T>;
	using pointer           = T*;
	using reference         = ref<T>;
	using iterator_category = std::random_access_iterator_tag;

	ptr() noexcept = default;
	explicit ptr(std::nullptr_t) noexcept {}
	template<class Other>
	constexpr explicit ptr(ptr<Other> const& /*other*/) noexcept {}
	// template<class Other,
	//          class = decltype(boost::multi::detail::implicit_cast<pointer>(std::declval<typename Other::pointer>()))>
	// constexpr explicit ptr(ptr<Other> const& /*other*/) noexcept {}  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions)  // NOSONAR(cpp:S1709)
	constexpr ptr(ptr const& /*other*/) = default;

	// vvv it is important that these two functions are device or device host functions
	// NOLINTNEXTLINE(fuchsia-overloaded-operator, fuchsia-trailing-return): this class simulates pointer
	constexpr auto operator*() const noexcept -> reference { return reference{}; }
	// NOLINTNEXTLINE(fuchsia-overloaded-operator, fuchsia-trailing-return): this class simulates pointer
	constexpr auto operator+(difference_type /*unused*/) const noexcept -> ptr { return *this; }
	// NOLINTNEXTLINE(fuchsia-overloaded-operator, fuchsia-trailing-return): this class simulates pointer

	auto operator+=(difference_type /*difference*/) noexcept -> ptr& { return *this; }
	// NOLINTNEXTLINE(fuchsia-overloaded-operator, fuchsia-trailing-return): this class simulates pointer
	auto operator++() noexcept -> ptr& { return operator+=(1); }
	// NOLINTNEXTLINE(fuchsia-overloaded-operator, fuchsia-trailing-return): this class simulates pointer
	friend auto operator-(ptr const& /*a*/, ptr const& /*b*/) noexcept -> difference_type { return 0; }
	// NOLINTNEXTLINE(fuchsia-overloaded-operator, fuchsia-trailing-return): this class simulates pointer
	auto operator==(ptr const& /*other*/) const noexcept -> bool { return true; }
	// NOLINTNEXTLINE(fuchsia-overloaded-operator, fuchsia-trailing-return): this class simulates pointer
	auto operator!=(ptr const& /*other*/) const noexcept -> bool { return false; }
	//  explicit operator T*() const{return &value;}
	// NOLINTNEXTLINE(fuchsia-overloaded-operator, fuchsia-trailing-return): this class simulates pointer
	auto operator->() const noexcept -> ptr const& { return *this; }
	// NOLINTNEXTLINE(fuchsia-trailing-return): this class simulates pointer
	//  friend auto to_address(ptr const& pointer) -> ptr {return pointer;}
	explicit operator bool() const noexcept { return false; }
	//  operator double*() const{return &value;}
	friend auto get_allocator(ptr const& /*self*/) noexcept { return std::allocator<value_type>{}; }
};

template<> double const ptr<double>::value       = 42.0;
template<> double const ptr<double const>::value = 42.0;

template<class T> class ref {
	friend class ptr<T>;
	friend class ref<T const>;
	ref() = default;

 public:
	//  explicit ref(ref<std::remove_const_t<T>> const& other) : p_{other.p_} {}
	// ~ref() = default;
	// auto operator=(ref const& other) -> ref& = delete;
	// constexpr ref(ref const& /*other*/) = delete;
	// constexpr ref(ref&& /*other*/) noexcept {}  // this is needed by nvcc, needs to be a device function for nvcc 11.2 and lower

	// auto operator=(ref     && other) noexcept -> ref& = delete;  // {*p_ = std::move(*other.p_); return *this;}

	constexpr operator T const&() const& { return ptr<T>::value; }  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions)  // NOSONAR(cpp:S1709)
	// NOLINTNEXTLINE(fuchsia-overloaded-operator): this class simulates a reference
	auto operator==(ref const& /*other*/) const { return true; }
	// NOLINTNEXTLINE(fuchsia-overloaded-operator): this class simulates a reference
	auto operator!=(ref const& /*other*/) const { return false; }
	using decay_t = std::decay_t<T>;
};

template<class T> struct allocator {
	using pointer    = ptr<T>;
	using value_type = T;
	auto allocate(std::size_t /*size*/) { return pointer{}; }
	void deallocate(pointer /*base*/, std::size_t /*size*/) {
		/*no-op;*/
	}
	//  std::true_type operator==(allocator const&){return {};}
	allocator() = default;
	template<class T2> explicit allocator(allocator<T2> const& /*other*/) {
		/*no-op;*/
	}
	template<class... Args>
	void construct(pointer /*location*/, Args const&... /*args*/) {
		/*no-op;*/
	}
	void destroy(pointer /*location*/) {
		/*no-op;*/
	}
};

// all these are optional, depending on the level of specialization needed
template<class Ptr, class T, class Size>
auto copy_n(Ptr /*first*/, Size /*count*/, ptr<T> result) {
	//  std::cerr<< "called Pointer-based copy_n(Ptr, n, fancy::ptr)" <<std::endl;
	return result;
}
template<class Ptr, class T, class Size>
auto copy_n(ptr<T> /*first*/, Size /*count*/, Ptr result) {
	//  std::cerr<< "called Pointer-based copy_n(fancy::ptr, n, Ptr)" <<std::endl;
	return result;
}
template<class T1, class T2, class Size>
auto copy_n(ptr<T1> /*first*/, Size /*count*/, ptr<T2> result) {
	//  std::cerr<< "called Pointer-based copy_n(fancy::ptr, n, fancy::ptr)" <<std::endl;
	return result;
}

}  // end namespace fancy

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  multi-fancy glue, where should this be?
//  In boost/multi/adaptors/MyFancyApaptor.hpp if anything, or in user code if it is very specialized
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

namespace boost::multi {

template<class It, class T>
auto copy(It first, It last, fancy::ptr<T> dest) {
	return copy(first, last, multi::array_iterator<T, 1, fancy::ptr<T>>{dest});
	//  std::cerr << "1D copy(it1D, it1D, it1D) with strides " << stride(first) << " " << stride(dest) << std::endl;
	//  return dest;
}

template<class It, class T>  // custom copy 1D (aka strided copy)
auto copy(It /*first*/, It /*last*/, multi::array_iterator<T, 1, fancy::ptr<T>> dest) {
	//  std::cerr << "1D copy(it1D, it1D, it1D) with strides " << stride(first) << " " << stride(dest) << std::endl;
	return dest;
}

template<class It, class T>  // custom copy 2D (aka double strided copy)
auto copy(It /*first*/, It /*last*/, multi::array_iterator<T, 2, fancy::ptr<T>> dest) {
	//  std::cerr<<"2D copy(It, It, it2D) with strides 1"<< first.stride() <<" "<< dest.stride() <<std::endl;
	return dest;
}

//  template<class Alloc, class It, class T> // custom copy 2D (aka double strided copy)
//  auto uninitialized_copy(Alloc&, It first, It last, multi::array_iterator<T, 2, fancy::ptr<T>> const& dest){
//     std::cerr << "2D uninitialized_copy(...) calls raw copy 2D" << std::endl;
//    return copy(first, last, dest);
//  }

}  // end namespace boost::multi

////////////////////////////////////////////////////////////////////////////////
// user code
////////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE(multi_fancy) {
	namespace multi = boost::multi;

	multi::array<double, 2, fancy::allocator<double>> arr({5, 5});
	BOOST_REQUIRE( arr.size() == 5 );
	BOOST_REQUIRE( arr[1][1] == arr[2][2] );

	multi::array<double, 2, fancy::allocator<double>> const arr2({0, 0});
	BOOST_REQUIRE( arr2.size() == 0 );
}
