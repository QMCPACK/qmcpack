// Copyright 2019-2023 Alfredo A. Correa

#include <boost/test/unit_test.hpp>

#include <multi/array.hpp>

#include <array>
#include <scoped_allocator>
#include <vector>

namespace multi = boost::multi;

template<class T = void>
class allocator1 {
	int* heap_ = nullptr;

	template<class> friend class allocator1;

 public:
	using value_type = T;

	allocator1() noexcept = delete;
	// NOLINTNEXTLINE(runtime/explicit)
	allocator1(int* heap) : heap_{heap} { assert(heap_); }  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions)
	template<class U> allocator1(allocator1<U> const& other) noexcept : heap_{other.heap_} {}  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions)

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
	friend auto operator==(allocator1 const& self, allocator1<U> const& other) noexcept { return self.heap_ == other.heap_; }
};

template<class T, class U>
auto operator!=(allocator1<T> const& self, allocator1<U> const& other) noexcept { return not(self == other); }

template<class T = void>
class allocator2 {
	std::int64_t* heap_ = nullptr;

	template<class> friend class allocator2;

 public:
	using value_type = T;

	allocator2() noexcept = delete;
	// NOLINTNEXTLINE(runtime/explicit)
	allocator2(std::int64_t* heap) : heap_{heap} { assert(heap_); }  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions)
	template<class U> allocator2(allocator2<U> const& other) noexcept : heap_{other.heap_} {}  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions)

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
	friend auto operator==(allocator2 const& self, allocator2<U> const& other) noexcept { return self.heap_ == other.heap_; }
};

template<class T, class U>
auto operator!=(allocator2<T> const& self, allocator2<U> const& other) noexcept { return not(self == other); }

#if not defined(__clang__) or not defined(__apple_build_version__)
BOOST_AUTO_TEST_CASE(scoped_allocator_vector) {
	std::int32_t heap1 = 0;
	std::int64_t heap2 = 0;

	{
		using InnerCont = std::vector<int, allocator2<int>>;
		using OuterCont = std::vector<InnerCont, std::scoped_allocator_adaptor<allocator1<InnerCont>, allocator2<int>>>;

		OuterCont cont({&heap1, &heap2});  // gives ambiguous construction in apple clang 14
		cont.resize(2);

		cont.resize(10);

		cont.back().resize(10);
		cont.back().resize(100);
		cont.back().resize(300);

		BOOST_TEST( heap1 == 1  );
		BOOST_TEST( heap2 == 1L );
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
		OuterCont cont({3, 4}, {&heap1, &heap2});  // gives ambiguous construction in apple clang 14

		cont[1][2].resize(10);
		cont[1][2].resize(100);
		cont[1][2].resize(200);

		BOOST_TEST( heap1 == 1  );
		BOOST_TEST( heap2 == 1L );
	}
}

BOOST_AUTO_TEST_CASE(scoped_allocator_array_vector_auto) {
	std::int32_t heap1 = 0;
	std::int64_t heap2 = 0;

	using InnerCont = std::vector<int, allocator2<int>>;
	using OuterCont = multi::array<InnerCont, 2, std::scoped_allocator_adaptor<allocator1<>, allocator2<>>>;

	{
		OuterCont cont({3, 4}, {&heap1, &heap2});  // gives ambiguous construction in apple clang 14

		cont[1][2].resize(10);
		cont[1][2].resize(100);
		cont[1][2].resize(200);

		BOOST_TEST( heap1 == 1  );
		BOOST_TEST( heap2 == 1L );
	}
}
#endif
