// Copyright 2020-2024 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_MULTI_ADAPTORS_FFTW_MEMORY_HPP
#define BOOST_MULTI_ADAPTORS_FFTW_MEMORY_HPP

#include <fftw3.h>

#include <cstddef>  // for std:::size_t

namespace boost::multi::fftw {

template<class T = void> struct allocator;

template<> struct allocator<void>{};

template<class T>
struct allocator {
	using value_type = T;
	using size_type  = std::size_t;

	static auto allocate(size_type n) -> T* {
		if(n == 0) {
			return nullptr;
		}
		void* ptr = fftw_malloc(sizeof(T) * n);
		if(ptr == nullptr) {
			throw std::bad_alloc{};
		}
		return static_cast<T*>(ptr);
	}
	static void deallocate(T* ptr, size_type n) {
		if(n != 0) {
			fftw_free(ptr);
		}
	}

	// constexpr auto operator==(allocator const& /*other*/) const -> bool { return true; }
	// constexpr auto operator!=(allocator const& /*other*/) const -> bool { return false; }
	constexpr friend auto operator==(allocator const& /*lhs*/, allocator const& /*rhs*/) { return true; }
	constexpr friend auto operator!=(allocator const& /*lhs*/, allocator const& /*rhs*/) { return false; }

	static constexpr auto max_size() noexcept { return std::numeric_limits<size_type>::max() / sizeof(T); }
};

template<class T, class U>
constexpr auto operator==(allocator<T> const& /*a*/, allocator<U> const& /*b*/) noexcept -> bool { return true; }

template<class T, class U>
constexpr auto operator!=(allocator<T> const& /*a*/, allocator<U> const& /*b*/) noexcept -> bool { return false; }

}  // namespace boost::multi::fftw
#endif  // BOOST_MULTI_ADAPTORS_FFTW_MEMORY_HPP
