// Copyright 2023-2024 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_MULTI_DETAIL_STATIC_ALLOCATOR_HPP
#define BOOST_MULTI_DETAIL_STATIC_ALLOCATOR_HPP

#include <boost/multi/detail/config/NODISCARD.hpp>
#include <boost/multi/detail/config/NO_UNIQUE_ADDRESS.hpp>

#include <array>
#include <cassert>
#include <cstddef>
#include <type_traits>

namespace boost::multi::detail {

template<class T, std::size_t N>
class static_allocator {  // NOSONAR(cpp:S4963) this allocator has special semantics
	bool dirty_ = false;

#ifdef _MSC_VER
	#pragma warning(push)
	#pragma warning(disable : 4324)  // Warning that the structure is padded due to the below
#endif

	BOOST_MULTI_NO_UNIQUE_ADDRESS alignas(T) std::array<std::byte, sizeof(T) * N> buffer_;

#ifdef _MSC_VER
	#pragma warning(pop)
#endif

 public:
	using value_type = T;
	using pointer    = T*;

	template<class TT> struct rebind {
		using other = static_allocator<TT, N>;
	};

	static constexpr auto max_size() noexcept -> std::size_t { return N; }

	static_allocator() = default;  // NOLINT(cppcoreguidelines-pro-type-member-init,hicpp-member-init) buffer_ is not initialized

	template<class TT, std::size_t NN>
	static_allocator(static_allocator<TT, NN> const& /*other*/) {  // NOLINT(hicpp-explicit-conversions,google-explicit-constructor) follow std::allocator
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

	static constexpr auto capacity() { return N; }

#ifdef _MSC_VER
#pragma warning( push )
#pragma warning( disable : 4068)  // bug in MSVC 14.2/14.3
#endif
	BOOST_MULTI_NODISCARD("because otherwise it will generate a memory leak")
	auto allocate([[maybe_unused]] std::size_t n) -> pointer {
		assert(n <= N);
		assert(!dirty_);  // do not attempt to resize a vector with static_allocator
		dirty_ = true;
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wcast-align"        // buffer_ is aligned as T
		return reinterpret_cast<pointer>(&buffer_);  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast)
#pragma GCC diagnostic pop
	}
#ifdef _MSC_VER
#pragma warning( pop ) 
#endif

	void deallocate(pointer /*ptr*/, [[maybe_unused]] std::size_t n) {
		assert(n <= N);
	}

	using is_always_equal = std::true_type;
};

template<class T, std::size_t N, class U>
constexpr auto operator==(static_allocator<T, N> const& /*a1*/, static_allocator<U, N> const& /*a2*/) noexcept { return true; }  // &a1 == &a2; }
// = delete;

template<class T, std::size_t N, class U>
auto operator!=(static_allocator<T, N> const& /*a1*/, static_allocator<U, N> const& /*a2*/) noexcept  // this is used *by the elements* when resizing a vector
{ return false; }                                                                                     // &a1 != &a2;}
// = delete

template<class T, std::size_t N, class U>
[[deprecated("don't swap dynamic container with static_allocator")]]
void swap(static_allocator<T, N>& a1, static_allocator<U, N>& a2) noexcept = delete;

}  // end namespace boost::multi::detail
#endif  // BOOST_MULTI_DETAIL_STATIC_ALLOCATOR_HPP
