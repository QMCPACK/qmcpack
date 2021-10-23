// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2021 Alfredo A. Correa

#ifndef MULTI_TUPLE_ZIP
#define MULTI_TUPLE_ZIP

#include<cassert>
#include<utility>

#include<tuple>

namespace boost {
namespace multi {
namespace detail {

//  Describe the type of a tuple with element I from each input tuple.
//  Needed to preserve the exact types from the input tuples.
//  template<std::size_t I, typename... Tuples>
//  using zip_tuple_at_index_t = std::tuple<std::tuple_element_t<I, std::decay_t<Tuples>>...>;

//  Collect all elements at index I from all input tuples as a new tuple.
//  template<std::size_t I, typename... Tuples>
//  zip_tuple_at_index_t<I, Tuples...> zip_tuple_at_index(Tuples && ...tuples) {
//    return {std::get<I>(std::forward<Tuples>(tuples))...};
//  }

template<class Tuple1, std::size_t... indices>
auto tuple_zip_impl(Tuple1&& t1, std::index_sequence<indices...> /*012*/) {
	return std::make_tuple(
		std::make_tuple(
			std::get<indices>(std::forward<Tuple1>(t1))
		)...
	);
}

template<class Tuple1, class Tuple2, std::size_t... indices>
auto tuple_zip_impl(Tuple1&& t1, Tuple2&& t2, std::index_sequence<indices...> /*012*/) {
	return std::make_tuple(
		std::make_tuple(
			std::get<indices>(std::forward<Tuple1>(t1)),
			std::get<indices>(std::forward<Tuple2>(t2))
		)...
	);
}

template<class Tuple1, class Tuple2, class Tuple3, std::size_t... indices>
auto tuple_zip_impl(Tuple1&& t1, Tuple2&& t2, Tuple3&& t3, std::index_sequence<indices...> /*012*/) {
	return	std::make_tuple(
		std::make_tuple(
			std::get<indices>(std::forward<Tuple1>(t1)),
			std::get<indices>(std::forward<Tuple2>(t2)),
			std::get<indices>(std::forward<Tuple3>(t3))
		)...
	);
}

template<class Tuple1, class Tuple2, class Tuple3, class Tuple4, std::size_t... indices>
auto tuple_zip_impl(Tuple1&& t1, Tuple2&& t2, Tuple3&& t3, Tuple4&& t4, std::index_sequence<indices...> /*012*/) {
	return std::make_tuple(
		std::make_tuple(
			std::get<indices>(std::forward<Tuple1>(t1)),
			std::get<indices>(std::forward<Tuple2>(t2)),
			std::get<indices>(std::forward<Tuple3>(t3)),
			std::get<indices>(std::forward<Tuple4>(t4))
		)...
	);
}

template<class Tuple1, class... Tuples>
auto tuple_zip(Tuple1&& t1, Tuples&&... ts) {
	return detail::tuple_zip_impl(
		std::forward<Tuple1>(t1), std::forward<Tuples>(ts)...,
		std::make_index_sequence<std::tuple_size<typename std::decay<Tuple1>::type>::value>()
	);
}

}  // end namespace detail
}  // end namespace multi
}  // end namespace boost

#endif


