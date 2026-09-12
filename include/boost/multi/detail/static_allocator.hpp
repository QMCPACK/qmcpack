// Copyright 2023-2026 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_MULTI_DETAIL_STATIC_ALLOCATOR_HPP
#define BOOST_MULTI_DETAIL_STATIC_ALLOCATOR_HPP
// #pragma once

#include "boost/multi/detail/config/NODISCARD.hpp"
#include "boost/multi/detail/config/NO_UNIQUE_ADDRESS.hpp"
#include "boost/multi/detail/implicit_cast.hpp"

#include <array>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <type_traits>

namespace boost::multi::detail {

// smallest unsigned integer type that can represent all values in [0, N]
template<std::uint64_t N>
using uint_for_bound_t =
	std::conditional_t<                             //
		(N <= 0xFFU), std::uint8_t,                 //
		std::conditional_t<                         //
			(N <= 0xFFFFU), std::uint16_t,          //
			std::conditional_t<                     //
				(N <= 0xFFFFFFFFU), std::uint32_t,  //
				std::uint64_t>>>;

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wpadded"
#endif

template<class T, typename Diff = typename std::pointer_traits<T*>::difference_type>
class offset_ptr {
	T*   ptr_;
	Diff offset_;

	template<class, typename> friend class offset_ptr;

 public:
	offset_ptr() = default;  // cppcheck-suppress uninitMemberVar
	constexpr explicit offset_ptr(T* ptr) : ptr_{ptr}, offset_{0} {}
	// cppcheck-suppress noExplicitConstructor ; nullptr should convert implicitly, like for raw pointers
	constexpr offset_ptr(std::nullptr_t) noexcept : ptr_{nullptr}, offset_{0} {}  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions,cppcoreguidelines-explicit-constructor,misc-explicit-constructor) match raw-pointer behavior

	offset_ptr(offset_ptr const&) = default;
	offset_ptr(offset_ptr&&)      = default;

	template<class U, std::enable_if_t<std::is_convertible_v<U*, T*>, int> = 0>  // NOLINT(modernize-use-constraints) for C++20
	// cppcheck-suppress noExplicitConstructor ; because underlying pointer is implicitly convertible
	constexpr /*mplct*/ offset_ptr(offset_ptr<U, Diff> const& other)  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions,cppcoreguidelines-explicit-constructor,misc-explicit-constructor) propagate implicitness of T*->U*
	: ptr_{other.ptr_}, offset_{other.offset_} {}

	template<class U, std::enable_if_t<std::is_constructible_v<T*, U*> && !std::is_convertible_v<U*, T*>, int> = 0>  // NOLINT(modernize-use-constraints) for C++20
	constexpr explicit offset_ptr(offset_ptr<U, Diff> const& other)
	: ptr_{static_cast<T*>(other.ptr_)}, offset_{other.offset_} {}

	auto operator=(offset_ptr const&) -> offset_ptr& = default;
	auto operator=(offset_ptr&&) -> offset_ptr&      = default;

	~offset_ptr() = default;

	using element_type      = typename std::pointer_traits<T*>::element_type;
	using difference_type   = Diff;
	using reference         = std::add_lvalue_reference_t<T>;  // T& for non-void T; void for T=void (T& would be ill-formed)
	using pointer           = T*;
	using value_type        = std::decay_t<T>;
	using iterator_category = std::random_access_iterator_tag;
#if defined(__cpp_lib_concepts) && (__cpp_lib_concepts >= 202002L)
	using iterator_concept = std::contiguous_iterator_tag;
#endif

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunknown-warning-option"
#pragma clang diagnostic ignored "-Wunsafe-buffer-usage"
#endif
	constexpr auto operator->() const noexcept -> pointer { return *this; }
	constexpr         operator T*() const noexcept {  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions,cppcoreguidelines-explicit-constructor,misc-explicit-constructor) implicit pointer-like conversion
        assert(ptr_ != nullptr || offset_ == 0);
        return ptr_ + offset_;  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
	}
#ifdef __clang__
#pragma clang diagnostic pop
#endif

	constexpr auto operator*() const noexcept -> reference { return *operator->(); }
	constexpr auto operator[](difference_type n) const -> reference { return *(*this + n); }

	constexpr auto operator-(offset_ptr const& other) const noexcept -> difference_type {
		assert(ptr_ == other.ptr_);
		return offset_ - other.offset_;
	}

	constexpr auto operator-(difference_type n) const noexcept -> offset_ptr { return offset_ptr{*this} -= n; }
	constexpr auto operator+(difference_type n) const noexcept -> offset_ptr { return offset_ptr{*this} += n; }

	friend constexpr auto operator+(difference_type n, offset_ptr const& ptr) noexcept -> offset_ptr { return ptr + n; }

	constexpr auto operator++() -> offset_ptr& {
		assert(ptr_ != nullptr);
		++offset_;
		return *this;
	}

	constexpr auto operator--() -> offset_ptr& {
		assert(ptr_ != nullptr);
		--offset_;
		return *this;
	}

	constexpr auto operator++(int) -> offset_ptr {
		auto ret = *this;
		++(*this);
		return ret;
	}

	constexpr auto operator--(int) -> offset_ptr {
		auto ret = *this;
		--(*this);
		return ret;
	}

	constexpr auto operator+=(difference_type n) noexcept -> offset_ptr& {
		assert(ptr_ != nullptr || n == 0);
		offset_ += n;
		return *this;
	}
	constexpr auto operator-=(difference_type n) noexcept -> offset_ptr& {
		return *this += -n;
	}

	friend constexpr auto operator==(offset_ptr const& lhs, offset_ptr const& rhs) noexcept -> bool {
		return lhs.ptr_ == rhs.ptr_ && lhs.offset_ == rhs.offset_;
	}

	friend constexpr auto operator!=(offset_ptr const& lhs, offset_ptr const& rhs) noexcept -> bool {
		return !(lhs == rhs);
	}

	friend constexpr auto operator<(offset_ptr const& lhs, offset_ptr const& rhs) noexcept -> bool {
		assert(lhs.ptr_ == rhs.ptr_);
		return lhs.offset_ < rhs.offset_;
	}

	friend constexpr auto operator<=(offset_ptr const& lhs, offset_ptr const& rhs) noexcept -> bool {
		assert(lhs.ptr_ == rhs.ptr_);
		return lhs.offset_ <= rhs.offset_;
	}

	friend constexpr auto operator>(offset_ptr const& lhs, offset_ptr const& rhs) noexcept -> bool {
		assert(lhs.ptr_ == rhs.ptr_);
		return lhs.offset_ > rhs.offset_;
	}

	friend constexpr auto operator>=(offset_ptr const& lhs, offset_ptr const& rhs) noexcept -> bool {
		assert(lhs.ptr_ == rhs.ptr_);
		return lhs.offset_ > rhs.offset_;
	}

	friend constexpr auto operator==(offset_ptr const& lhs, std::nullptr_t const& rhs) noexcept -> bool {
		return lhs.ptr_ == rhs;
	}
	friend constexpr auto operator!=(offset_ptr const& lhs, std::nullptr_t const& rhs) noexcept -> bool {
		return !(lhs == rhs);
	}
};

template<class T, std::size_t N, typename Ptr = offset_ptr<T>>
class static_allocator {  // NOSONAR(cpp:S4963) this allocator has special semantics
#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable : 4324)  // Warning that the structure is padded due to the below
#endif

	// #if defined(__clang__)
	// #pragma clang diagnostic push
	// #pragma clang diagnostic ignored "-Wpadded"
	// #endif

	BOOST_MULTI_NO_UNIQUE_ADDRESS alignas(T) std::array<std::byte, sizeof(T) * N> buffer_;

	// #if defined(__clang__)
	// #pragma clang diagnostic pop
	// #endif

#ifdef _MSC_VER
#pragma warning(pop)
#endif

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable : 4820)  // warning C4820: 'boost::multi::detail::static_allocator<main::T,32>': '3' bytes padding added after data member 'boost::multi::detail::static_allocator<main::T,32>::dirty_' [C:\Gitlab-Runner\builds\t3_1sV2uA\0\correaa\boost-multi\build\test\allocator.cpp.x.vcxproj]
#endif
	bool dirty_ = false;
#ifdef _MSC_VER
#pragma warning(pop)
#endif

 public:
	using value_type = T;
	using pointer    = Ptr;

	template<class TT> struct rebind {
		using other = static_allocator<TT, N>;
	};

	// using size_type = std::conditional_t<(N <= 0xFFFF), std::uint16_t, std::uint32_t>;

	using size_type = std::conditional_t<(N <= 0xFFFF), std::uint16_t, std::uint32_t>;

	static constexpr auto max_size() noexcept -> size_type { return N; }

	static_allocator() = default;  // NOLINT(cppcoreguidelines-pro-type-member-init,hicpp-member-init) buffer_ is not initialized

	template<class TT, std::size_t NN>
	explicit static_allocator(static_allocator<TT, NN> const& /*other*/) {  // NOLINT(hicpp-explicit-conversions,google-explicit-constructor) follow std::allocator  // NOSONAR
		// static_assert(sizeof(T) == sizeof(TT));
		static_assert(NN == N);
	}

	static_allocator(static_allocator const& /*other*/)  // std::vector makes a copy right away
	// = default;  // this copies the internal buffer
	{}

	// [[deprecated("don't move dynamic container with static_allocator")]]
	static_allocator(static_allocator&& /*other*/)  // this is called *by the elements* during move construction of a vector
													// = delete;
													// {throw std::runtime_error("don't move dynamic container with static_allocator");}  // this is called *by the elements* during move construction of a vector
		noexcept {}
	// noexcept {std::memmove(buffer_.data(), other.buffer_.data(), sizeof(T)*N);}
	// noexcept : buffer_{std::move(other.buffer_)} {}
	// noexcept = default;

	[[deprecated("don't move dynamic container with static_allocator")]]
	auto operator=(static_allocator const& /*other*/) -> static_allocator& = delete;

	[[deprecated("don't move dynamic container with static_allocator")]] auto operator=(static_allocator&& other) -> static_allocator& = delete;

	~static_allocator() = default;

	auto select_on_container_copy_construction() noexcept -> static_allocator = delete;
	// {return static_allocator{};}

	using propagate_on_container_move_assignment = std::false_type;  // this forces to call move assignment of the allocator by std::vector
	using propagate_on_container_copy_assignment = std::false_type;
	using propagate_on_container_swap            = std::false_type;

	static constexpr auto capacity() -> size_type { return N; }

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable : 4068)  // bug in MSVC 14.2/14.3
#endif
	BOOST_MULTI_NODISCARD("because otherwise it will generate a memory leak")
	auto allocate([[maybe_unused]] size_type n) -> pointer {
		assert(n <= N);
		assert(!dirty_);  // do not attempt to resize a vector with static_allocator
						  // dirty_ = true;
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wcast-align"            // buffer_ is aligned as T
		return pointer(reinterpret_cast<T*>(&buffer_));  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast)
#pragma GCC diagnostic pop
	}
#ifdef _MSC_VER
#pragma warning(pop)
#endif

	static void deallocate(pointer /*ptr*/, [[maybe_unused]] size_type n) {
		assert(n <= N);
	}

	using is_always_equal = std::true_type;
};

#ifdef __clang__
#pragma clang diagnostic pop
#endif

template<class T, std::size_t N, class U>
constexpr auto operator==(static_allocator<T, N> const& /*a1*/, static_allocator<U, N> const& /*a2*/) noexcept { return true; }  // &a1 == &a2; }
// = delete;

template<class T, std::size_t N, class U>
auto operator!=(static_allocator<T, N> const& /*a1*/, static_allocator<U, N> const& /*a2*/) noexcept  // this is used *by the elements* when resizing a vector
{ return false; }                                                                                     // &a1 != &a2;}
// = delete

template<class T, std::size_t N, class U>
[[deprecated("don't swap dynamic container with static_allocator")]]
void swap(static_allocator<T, N>& left, static_allocator<U, N>& right) noexcept = delete;

}  // end namespace boost::multi::detail
#endif  // BOOST_MULTI_DETAIL_STATIC_ALLOCATOR_HPP
