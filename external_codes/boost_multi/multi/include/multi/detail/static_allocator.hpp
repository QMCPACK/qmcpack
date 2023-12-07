// Copyright 2023 Alfredo A. Correa

#ifndef MULTI_DETAIL_STATIC_ALLOCATOR_HPP
#define MULTI_DETAIL_STATIC_ALLOCATOR_HPP

#include "../config/NODISCARD.hpp"
#include "../config/NO_UNIQUE_ADDRESS.hpp"

#include <cassert>

namespace boost::multi::detail {

template<class T, std::size_t N>
class static_allocator {
	bool dirty_ = false;
	MULTI_NO_UNIQUE_ADDRESS alignas(T) std::array<std::byte, sizeof(T) * N> buffer_;

 public:
	using value_type = T;
	using pointer    = T*;

	template<class TT> struct rebind {
		using other = static_allocator<TT, N>;
	};

	static constexpr auto max_size() noexcept -> std::size_t { return N; };

	static_allocator() = default;  // NOLINT(cppcoreguidelines-pro-type-member-init,hicpp-member-init) buffer_ is not initialized

	static_allocator(static_allocator const& /*other*/) // std::vector makes a copy right away
	// = default;  // this copies the internal buffer
	{}
	[[deprecated("don't move dynamic container with static_allocator")]]
	static_allocator(static_allocator&& /*other*/)  // this is called *by the elements* during move construction of a vector
	// = delete;
	// {throw std::runtime_error("don't move dynamic container with static_allocator");}  // this is called *by the elements* during move construction of a vector
	noexcept {}
	// noexcept {std::memmove(buffer_.data(), other.buffer_.data(), sizeof(T)*N);}
	// noexcept : buffer_{std::move(other.buffer_)} {}
	// noexcept = default;

	[[deprecated("don't move dynamic container with static_allocator")]]
	auto operator=(static_allocator const& /*other*/) -> static_allocator&
	= delete;

	[[deprecated("don't move dynamic container with static_allocator")]] auto operator=(static_allocator&& other) -> static_allocator&
	= delete;

	~static_allocator() = default;

	auto select_on_container_copy_construction() noexcept -> static_allocator
	= delete;
	// {return static_allocator{};}

	using propagate_on_container_move_assignment = std::false_type;  // this forces to call move assignment of the allocator by std::vector
	using propagate_on_container_copy_assignment = std::false_type;
	using propagate_on_container_swap            = std::false_type;

	static constexpr auto capacity() { return N; }

	NODISCARD("because otherwise it will generate a memory leak")
	auto allocate([[maybe_unused]] std::size_t n) -> pointer {
		assert(n <= N);
		assert(not dirty_);  // do not attempt to resize a vector with static_allocator
		dirty_ = true;
		return reinterpret_cast<pointer>(buffer_.data());  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast)
	}
	void deallocate(pointer /*ptr*/, [[maybe_unused]] std::size_t n) {
		assert(n <= N);
	}

	using is_always_equal = std::true_type;
};

template<class T, std::size_t N, class U>
constexpr auto operator==(static_allocator<T, N> const& /*a1*/, static_allocator<U, N> const& /*a2*/) noexcept
{ return true; } // &a1 == &a2; }
// = delete;

template <class T, std::size_t N, class U>
auto operator!=(static_allocator<T, N> const& /*a1*/, static_allocator<U, N> const& /*a2*/) noexcept // this is used *by the elements* when resizing a vector
{ return false; } // &a1 != &a2;}
// = delete

template <class T, std::size_t N, class U>
[[deprecated("don't swap dynamic container with static_allocator")]] void swap(static_allocator<T, N>& a1, static_allocator<U, N>& a2) = delete;

}  // end namespace boost::multi::detail
#endif  // MULTI_DETAIL_STATIC_ALLOCATOR_HPP
