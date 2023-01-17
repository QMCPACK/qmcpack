// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2018-2022 Alfredo A. Correa

#ifndef MPI3_DETAIL_TUPLE_OFFSET_HPP
#define MPI3_DETAIL_TUPLE_OFFSET_HPP

#include <tuple>

namespace boost::mpi3::detail {

template<std::size_t I, typename Tuple>
[[deprecated]] constexpr std::size_t tuple_offset_aux() {  // from https://stackoverflow.com/questions/70647441/how-to-determine-the-offset-of-an-element-of-a-tuple-at-compile-time
	static_assert(not std::is_reference_v<std::tuple_element_t<I, Tuple>>);
	union u {
		constexpr u() : a{} {}  // GCC bug needs a constructor definition
		char  a[sizeof(Tuple)];  // NOLINT(cppcoreguidelines-avoid-c-arrays,cppcoreguidelines-pro-bounds-pointer-arithmetic,hicpp-avoid-c-arrays,modernize-avoid-c-arrays)  pretty low level stuff
		Tuple t;
	} x;
	auto off =
		std::find_if(
			static_cast<char*>(x.a), static_cast<char*>(x.a) + sizeof(Tuple),  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic) pretty low level stuff
			[&x](char const& e) { std::addressof(e) == std::addressof(std::get<I>(x.t)); }
		) -
		static_cast<char*>(x.a);
	// std::size_t off = 0;
	// while(static_cast<void*>(x.a + off) != std::addressof(std::get<I>(x.t))) {++off;}  // NOLINT(cppcoreguidelines-pro-type-union-access,cppcoreguidelines-pro-bounds-pointer-arithmetic)
	return off;
}

template<std::size_t I, typename Tuple>
std::size_t element_offset() {
	static std::size_t const ret = [] {
		static_assert(!std::is_reference_v<std::tuple_element_t<I, Tuple>>);
		union u {
			constexpr u() : a{} {}  // GCC bug needs a constructor definition
			char  a[sizeof(Tuple)];  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays)
			Tuple t;
		} x;
		return reinterpret_cast<char*>(std::addressof(std::get<I>(x.t))) - x.a;  // NOLINT(cppcoreguidelines-pro-type-union-access,cppcoreguidelines-pro-type-reinterpret-cast)
	}();
	return ret;
}

template<std::size_t Iter, typename Tuple>
struct tuple_offset : std::integral_constant<std::size_t, tuple_offset_aux<Iter, Tuple>()> {};

template<std::size_t Iter, typename Tuple>
constexpr std::size_t tuple_offset_v = tuple_offset<Iter, Tuple>::value;

}  // end namespace boost::mpi3::detail
#endif
