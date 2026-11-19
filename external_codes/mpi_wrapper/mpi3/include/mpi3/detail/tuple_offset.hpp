// Copyright 2018-2025 Alfredo A. Correa

#ifndef MPI3_DETAIL_TUPLE_OFFSET_HPP
#define MPI3_DETAIL_TUPLE_OFFSET_HPP

#include <iterator>
#include <tuple>

namespace boost::mpi3::detail {

template <std::size_t I, typename Tuple>
constexpr auto element_offset() -> std::size_t {
    static_assert(!std::is_reference_v<std::tuple_element_t<I, Tuple> >);
    union u {
        constexpr u() : a{} {}
		char a[sizeof(Tuple)];  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) can't use std::array here in nvcc 24
        Tuple t;
    } x;
	using std::get;
	static_assert(std::is_reference_v<decltype(get<I>(x.t))>);  // NOLINT(cppcoreguidelines-pro-type-union-access)
	auto const* const p = std::addressof(get<I>(x.t));  // NOLINT(cppcoreguidelines-pro-type-union-access)
	std::size_t ret = 0;
    while(static_cast<void*>(x.a + ret) != p) {  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic,cppcoreguidelines-pro-type-union-access,altera-unroll-loops,altera-id-dependent-backward-branch)
		++ret;
	}
	return ret;
}

template<std::size_t I, typename Tuple>
struct tuple_offset : std::integral_constant<std::size_t, element_offset<I, Tuple>()> {};

template<std::size_t I, typename Tuple>
constexpr std::size_t tuple_offset_v = tuple_offset<I, Tuple>::value;

}  // end namespace boost::mpi3::detail
#endif
