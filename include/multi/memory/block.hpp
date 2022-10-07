// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2018-2022 Alfredo A. Correa

#ifndef BOOST_MULTI_MEMORY_BLOCK_HPP
#define BOOST_MULTI_MEMORY_BLOCK_HPP

#include <cassert>   // for assert
#include <cstddef>   // for nullptr_t, size_t
#include <iterator>  // for distance
#include <memory>    // for pointer_traits

#ifdef USE_BOOST_MULTI_MEMORY_MALLOC
#include <malloc.h>
#endif

namespace boost {
namespace multi {
namespace memory {

#if(__cpp_lib_byte<201603)
enum class byte : unsigned char {};
#else
using byte = std::byte;
#endif

template<
	typename Ptr = byte*,
	typename Diff = typename std::pointer_traits<Ptr>::difference_type
> 
struct block;

namespace detail {
template<typename Ptr, typename Difference>
struct basic_block {
	using pointer = typename std::pointer_traits<Ptr>::pointer; // basic_block? block<Ptr>?
	using element_type = typename std::pointer_traits<Ptr>::element_type;
	using difference_type = Difference;//typename std::pointer_traits<Ptr>::difference_type;
	using size_type = difference_type;

 public:
	typename basic_block::pointer start_;
	size_type size_;
	constexpr basic_block() = default;
	constexpr basic_block(pointer p, size_type size) : start_{p}, size_{size}{}
	constexpr basic_block(std::nullptr_t p, size_type s): start_{p}, size_{s}{assert(!s);}
#ifdef USE_BOOST_MULTI_MEMORY_MALLOC
	template<class OtherPtr>
	constexpr basic_block(OtherPtr p)
	: start_{p}, size_{static_cast<size_type>(malloc_usable_size(p))}{assert(size_);}
#endif
	constexpr size_type size() const{return size_;}
	constexpr operator Ptr const&() const{return start_;}
	constexpr bool contains(typename basic_block::pointer p) const{
		using std::distance;
		difference_type d = distance(start_, p);
		return d >= 0 and d < size_;
	}
	friend constexpr size_type size(basic_block const& self){return self.size();}
};
}

template<typename Ptr, typename Diff>
struct block : detail::basic_block<Ptr, Diff> {
	using detail::basic_block<Ptr, Diff>::basic_block;
};

template<typename T, typename Diff>
struct block<T*, Diff> : detail::basic_block<T*, Diff> {
	using detail::basic_block<T*, Diff>::basic_block;
	template<std::size_t N>
	constexpr explicit block(T(&t)[N]) : detail::basic_block<T*, Diff>{t, N}{}
};

#if (__cpp_deduction_guides >= 201703)
template<class Ptr, class Size> block(Ptr p, Size s) -> block<Ptr, Size>;
#endif

#ifdef USE_BOOST_MULTI_MEMORY_MALLOC
#if (__cpp_deduction_guides >= 201703)
template<class Ptr> 
block(Ptr p) -> block<Ptr, typename std::pointer_traits<Ptr>::difference_type>;
#endif
template<typename Ptr>
constexpr size_t size(Ptr const& p){return block<Ptr>(p).size();}//size();}//static_cast<block<Ptr> const&>(p).size();}
#endif

}}}

#if _TEST_BOOST_MULTI_MEMORY_BLOCK

#include<cassert>
#include<vector>

namespace multi = boost::multi;
namespace memory = multi::memory;

int main() {

	char A[1024];  // flawfinder: ignore testing legacy type
	memory::block<char*> a{A};

	std::vector<char> V(1000);
	memory::block<char*> v(V.data(), 500);
	memory::block<char*> v2 = v; (void)v2;

	memory::byte B[1024];
	memory::block<> b{B};

	{
		std::vector<char> V(1000);
		memory::block<char*, std::integral_constant<std::size_t, 500>> v(V.data(), {});
		memory::block<char*> v2 = v; (void)v2;
	}
	#if (__cpp_deduction_guides >= 201703)
	{
		std::vector<char> V(1000);
		memory::block v(V.data(), std::integral_constant<std::size_t, 500>{});
		memory::block<char*, std::ptrdiff_t> v2 = v; (void)v2;
		assert( v == v2 );
	}
	#endif
#ifdef USE_BOOST_MULTI_MEMORY_MALLOC
	{
		char* p = new char[200];
		memory::block<char*> bp = p; // non standard
		assert( bp.size() >= 200 );
		using multi::memory::size;
		assert( size(p) >= 200 );
	}
	{
		memory::block<char*> b = new char[200];
		*b = 'a';
		*(b + 1) = 'b';
		assert( *b == 'a' );
		assert( b[0] == 'a' );
		assert( b[1] == 'b' );
		assert( size(b) >= 200 );
		delete b;
	}
	#if (__cpp_deduction_guides >= 201703)
	{
		memory::block b = new char[200];
		assert( size(b) >= 200 );
		delete b;
	}
	#endif
#endif

}
#endif
#endif
