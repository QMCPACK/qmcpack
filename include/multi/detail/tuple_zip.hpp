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

template<class Tuple1, std::size_t... Indices>
auto tuple_zip_impl(Tuple1&& t1, std::index_sequence<Indices...> /*012*/) {
	return std::make_tuple(
		std::make_tuple(
			std::get<Indices>(std::forward<Tuple1>(t1))
		)...
	);
}

template<class Tuple1, class Tuple2, std::size_t... Is>
auto tuple_zip_impl(Tuple1&& t1, Tuple2&& t2, std::index_sequence<Is...> /*012*/) {
	return std::make_tuple(
		std::make_tuple(
			std::get<Is>(std::forward<Tuple1>(t1)),
			std::get<Is>(std::forward<Tuple2>(t2))
		)...
	);
}

template<class Tuple1, class Tuple2, class Tuple3, std::size_t... Is>
auto tuple_zip_impl(Tuple1&& t1, Tuple2&& t2, Tuple3&& t3, std::index_sequence<Is...> /*012*/) {
	return	std::make_tuple(
		std::make_tuple(
			std::get<Is>(std::forward<Tuple1>(t1)),
			std::get<Is>(std::forward<Tuple2>(t2)),
			std::get<Is>(std::forward<Tuple3>(t3))
		)...
	);
}

template<class Tuple1, class Tuple2, class Tuple3, class Tuple4, std::size_t... Is>
auto tuple_zip_impl(Tuple1&& t1, Tuple2&& t2, Tuple3&& t3, Tuple4&& t4, std::index_sequence<Is...> /*012*/) {
	return std::make_tuple(
		std::make_tuple(
			std::get<Is>(std::forward<Tuple1>(t1)),
			std::get<Is>(std::forward<Tuple2>(t2)),
			std::get<Is>(std::forward<Tuple3>(t3)),
			std::get<Is>(std::forward<Tuple4>(t4))
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
