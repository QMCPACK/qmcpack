// Copyright 2019-2024 Alfredo A. Correa
// Copyright 2024 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/multi/array.hpp>

#include <array>
#include <scoped_allocator>
#include <vector>

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

template<class T = void>
class allocator1 {
	int* heap_ = nullptr;

	template<class> friend class allocator1;

 public:
	using value_type = T;

	allocator1() noexcept = delete;
	// NOLINTNEXTLINE(runtime/explicit)
	allocator1(int* heap) : heap_{heap} { assert(heap_); }  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions)  // NOSONAR(cpp:S1709) mimic memory resource syntax (pass pointer)
	template<class U> allocator1(allocator1<U> const& other) noexcept : heap_{other.heap_} {}  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions)  // NOSONAR(cpp:S1709) allocator conversions are not explicit

	auto allocate(std::size_t n) {
		if(n == 0) {
			return static_cast<value_type*>(nullptr);
		}
		if(heap_ == nullptr) {
			throw std::bad_alloc{};
		}  // this cuts branches with UB (null deref) for the sanitizer
		++*heap_;
		return static_cast<value_type*>(::operator new(n * sizeof(value_type)));
	}
	void deallocate(value_type* ptr, std::size_t n) noexcept {
		if(n == 0) {
			return;
		}
		--*heap_;
		::operator delete(ptr);
	}
	template<class U>
	friend auto operator==(allocator1 const& self, allocator1<U> const& other) noexcept -> bool { return self.heap_ == other.heap_; }
	template<class U>
	friend auto operator!=(allocator1 const& self, allocator1<U> const& other) noexcept -> bool { return self.heap_ != other.heap_; }
};

template<class T, class U>
auto operator!=(allocator1<T> const& self, allocator1<U> const& other) noexcept -> bool { return ! (self == other); }

template<class T, class U>
auto operator==(allocator1<T> const& self, allocator1<U> const& other) noexcept -> bool { return   (self == other); }

template<class T = void>
class allocator2 {
	std::int64_t* heap_ = nullptr;

	template<class> friend class allocator2;

 public:
	using value_type = T;

	allocator2() noexcept = default;
	// NOLINTNEXTLINE(runtime/explicit)
	allocator2(std::int64_t* heap) : heap_{heap} { assert(heap_); }  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions)  // NOSONAR(cpp:S1709) mimic memory resource syntax (pass pointer)
	template<class U> allocator2(allocator2<U> const& other) noexcept : heap_{other.heap_} {}  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions)  // NOSONAR(cpp:S1709) allocator conversions are not explicit

	auto allocate(std::size_t n) {
		if(n == 0) {
			return static_cast<value_type*>(nullptr);
		}
		if(heap_ == nullptr) {
			throw std::bad_alloc{};
		}  // this cuts branches with UB (null deref) for the sanitizer
		++*heap_;
		return static_cast<value_type*>(::operator new(n * sizeof(value_type)));
	}

	void deallocate(value_type* ptr, std::size_t n) noexcept {
		if(n == 0) {
			return;
		}
		--*heap_;
		::operator delete(ptr);
	}

	template<class U>
	friend auto operator==(allocator2 const& self, allocator2<U> const& other) noexcept -> bool { return self.heap_ == other.heap_; }
	template<class U>
	friend auto operator!=(allocator2 const& self, allocator2<U> const& other) noexcept -> bool { return self.heap_ != other.heap_; }
};

template<class T, class U>
auto operator!=(allocator2<T> const& self, allocator2<U> const& other) noexcept -> bool {
	return ! (self == other);
}

template<class T, class U>
auto operator==(allocator2<T> const& self, allocator2<U> const& other) noexcept -> bool {
	return (self == other);
}

BOOST_AUTO_TEST_CASE(scoped_allocator_vector) {
	std::int32_t heap1 = 0;
	std::int64_t heap2 = 0;

	{
		using InnerCont = std::vector<int, allocator2<int>>;
		using OuterCont =
			std::vector<
				InnerCont,
				std::scoped_allocator_adaptor<
					allocator1<InnerCont>,
					allocator2<int>
				>
			>
		;

		// OuterCont cont({&heap1, &heap2});  // gives ambiguous construction in libc++
		OuterCont cont({&heap1, allocator2<int>{&heap2}});

		cont.resize(2);

		cont.resize(10);

		cont.back().resize(10);
		cont.back().resize(100);
		cont.back().resize(300);

	// these values are depdenent on the implementation of std::vector
	#if !defined(_MSC_VER)
		BOOST_TEST( heap1 == 1  );
		BOOST_TEST( heap2 == 1L );
	#endif
	}

	BOOST_TEST( heap1 == 0 );
	BOOST_TEST( heap2 == 0 );
}

BOOST_AUTO_TEST_CASE(scoped_allocator_array_vector) {
	std::int32_t heap1 = 0;
	std::int64_t heap2 = 0;

	using InnerCont = std::vector<int, allocator2<int>>;
	using OuterCont = multi::array<InnerCont, 2, std::scoped_allocator_adaptor<allocator1<InnerCont>, allocator2<int>>>;

	{
		OuterCont cont(
		#ifdef _MSC_VER  // problem with MSVC 14.3 c++17
			multi::extensions_t<2>
		#endif
			{3, 4},
			{&heap1, allocator2<int>{&heap2}}  // without allocator2<>{...} gives ambiguous construction in libc++
		);

		cont[1][2].resize(10);
		cont[1][2].resize(100);
		cont[1][2].resize(200);

	// these values are depdenent on the implementation of std::vector
	#if !defined(_MSC_VER)
		BOOST_TEST( heap1 == 1  );
		BOOST_TEST( heap2 == 1L );
	#endif
	}
}

// vvv this cases confuse gcc (and MSVC?)
// BOOST_AUTO_TEST_CASE(scoped_allocator_array_vector_auto) {
//  std::int32_t heap1 = 0;
//  std::int64_t heap2 = 0;

//  using InnerCont = std::vector<int, allocator2<int>>;
//  using OuterCont = multi::array<InnerCont, 2, std::scoped_allocator_adaptor<allocator1<>, allocator2<>>>;

//  {
//    OuterCont cont({3, 4}, {&heap1, allocator2<>{&heap2}});  // without allocator2<>{...} gives ambiguous construction in libc++

//    cont[1][2].resize( 10);
//    cont[1][2].resize(100);
//    cont[1][2].resize(200);

//    BOOST_TEST( heap1 == 1  );
//  // these values are depdenent on the implementation of std::vector
//  #if !defined(_MSC_VER)
//    BOOST_TEST( heap2 ==  1L );
//  #endif
//  }
// }

// BOOST_AUTO_TEST_CASE(scoped_allocator_array_array_auto) {
//  std::int32_t heap1 = 0;
//  std::int64_t heap2 = 0;

//  using InnerCont = multi::array<int, 2, allocator2<int>>;
//  using OuterCont = multi::array<InnerCont, 2, std::scoped_allocator_adaptor<allocator1<>, allocator2<>>>;

//  {
//    OuterCont cont({3, 4}, {&heap1, allocator2<>{&heap2}});  // without allocator2<>{...} gives ambiguous construction in libc++

//    cont[1][2].reextent({ 10,  10});
//    cont[1][2].reextent({100, 100});
//    cont[1][2].reextent({200, 200});

//    BOOST_TEST( heap1 == 1  );
//    BOOST_TEST( heap2 == 1L );
//  }
// }
