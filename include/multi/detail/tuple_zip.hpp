// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2021-2022 Alfredo A. Correa

#ifndef MULTI_DETAIL_TUPLE_ZIP_HPP
#define MULTI_DETAIL_TUPLE_ZIP_HPP

#include<cassert>
#include<utility>

#include<tuple>  // for deprecated functions

namespace boost::multi {  // NOLINT(modernize-concat-nested-namespaces) keep c++14 compat
namespace detail {

template<class... Ts> class tuple;

template<class T0, class... Ts> class tuple<T0, Ts...> {  // NOLINT(cppcoreguidelines-special-member-functions,hicpp-special-member-functions)
	T0 head_;
	tuple<Ts...> tail_;

 public:
	constexpr auto head() const& -> T0 const& {return head_;}
	constexpr auto tail() const& -> tuple<Ts...> const& {return tail_;}

	constexpr auto head() && -> T0&& {return std::move(head_);}
	constexpr auto tail() && -> tuple<Ts...>&& {return std::move(tail_);}

	constexpr auto head() & -> T0& {return head_;}
	constexpr auto tail() & -> tuple<Ts...>& {return tail_;}

	tuple() = default;
	tuple(tuple const&) = default;

	template<class TT0, class... TTs>
	[[deprecated]]
	// cppcheck-suppress noExplicitConstructor ; for compatibility with QMC to be deprecated
	constexpr tuple(std::tuple<TT0, TTs...> const& other)  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions)
	: tuple{std::apply([](auto const&... es) noexcept {return tuple(es...);}, other)}  // paren to allow narrowing conversions
	{}

	constexpr auto operator=(tuple const& other) -> tuple& {
		if(this == &other) {return *this;}
		head_ = other.head_;
		tail_ = other.tail_;
		return *this;
	}

	template<class TT0>
	// cppcheck-suppress noExplicitConstructor ; for compatibility with QMC to be deprecated
	[[deprecated]] constexpr tuple(TT0 t0)  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions)
	: head_{std::move(t0)}, tail_{} {}

	constexpr explicit tuple(T0 t0, Ts... ts)
	: head_{std::move(t0)}, tail_{std::move(ts)...} {}

	// cppcheck-suppress noExplicitConstructor ; allow bracket init in function argument // NOLINTNEXTLINE(runtime/explicit)
	constexpr tuple(T0 t0, tuple<Ts...> sub) : head_{std::move(t0)}, tail_{std::move(sub)} {}

	constexpr auto operator==(tuple const& other) const -> bool {
		return head_ == other.head_ and tail_ == other.tail_;
	}
	constexpr auto operator!=(tuple const& other) const -> bool {
		return head_ != other.head_ or  tail_ != other.tail_;
	}

	constexpr auto operator< (tuple const& other) const {
		if(head_ < other.head_) {return true;}
		if(other.head_ < head_) {return false;}
		return tail_ < other.tail_;
	}
	constexpr auto operator> (tuple const& other) const {
		if(head_ > other.head_) {return true;}
		if(other.head_ > head_) {return false;}
		return tail_ > other.tail_;
	}
};

template<> class tuple<> {
 public:
	constexpr auto operator==(tuple const& /*other*/) const -> bool {return true ;}
	constexpr auto operator!=(tuple const& /*other*/) const -> bool {return false;}

	constexpr auto operator< (tuple const& /*other*/) const {return false;}
	constexpr auto operator> (tuple const& /*other*/) const {return false;}
};

template<class T0, class... Ts> tuple(T0, tuple<Ts...>) -> tuple<T0, Ts...>;

template<class T0, class... Ts> constexpr auto mk_tuple(T0 t0, Ts... ts) {
	return tuple<T0, Ts...>(std::move(t0), std::move(ts)...);
}

template<class T0, class Tuple> struct tuple_prepend;

template<class T0, class... Ts>
struct tuple_prepend<T0, tuple<Ts...>> {
	using type = tuple<T0, Ts...>;
};

template<class T0, class Tuple>
using tuple_prepend_t = typename tuple_prepend<T0, Tuple>::type;

template<class T0, class... Ts>
constexpr auto head(tuple<T0, Ts...> const& t) ->  decltype(auto) {
	return t.head();
}

template<class T0, class... Ts>
constexpr auto head(tuple<T0, Ts...> && t) -> decltype(auto) {
	return std::move(t.head());
}

template<class T0, class... Ts>
constexpr auto head(tuple<T0, Ts...> & t) -> decltype(auto) {
	return t.head();
}

template<class T0, class... Ts>
constexpr auto tail(tuple<T0, Ts...> const& t) -> decltype(auto) {return t.tail();}

template<class T0, class... Ts>
constexpr auto tail(tuple<T0, Ts...>     && t) -> decltype(auto) {return std::move(t).tail();}

template<class T0, class... Ts>
constexpr auto tail(tuple<T0, Ts...>      & t) -> decltype(auto) {return t.tail();}

#if defined __NVCC__
    #ifdef __NVCC_DIAG_PRAGMA_SUPPORT__
        #pragma nv_diagnostic push
        #pragma nv_diag_suppress = implicit_return_from_non_void_function
    #else
        #pragma    diagnostic push
        #pragma    diag_suppress = implicit_return_from_non_void_function
    #endif
#elif defined __NVCOMPILER
    #pragma    diagnostic push
    #pragma    diag_suppress = implicit_return_from_non_void_function
#endif
template<std::size_t N, class T0, class... Ts>
constexpr auto get(tuple<T0, Ts...> const& t) -> auto const& {
	if constexpr(N == 0) {
		return t.head();
	} else {
		return get<N-1>(t.tail());
	}
}

template<std::size_t N, class T0, class... Ts>
constexpr auto get(tuple<T0, Ts...>& t) -> auto& {
	if constexpr(N == 0) {
		return t.head();
	} else {
		return get<N-1>(t.tail());
	}
}

template<std::size_t N, class T0, class... Ts>
constexpr auto get(tuple<T0, Ts...>&& t) -> auto&& {
	if constexpr(N == 0) {
		return std::move(t).head();
	} else {
		return get<N-1>(std::move(t.tail()));
	}
}
#if defined __NVCC__
    #ifdef __NVCC_DIAG_PRAGMA_SUPPORT__
        #pragma nv_diagnostic pop
    #else
        #pragma    diagnostic pop
    #endif
#elif defined __NVCOMPILER
    #pragma    diagnostic pop
#endif

}  // end namespace detail
}  // end namespace boost::multi

namespace std {  // NOLINT(cert-dcl58-cpp) define stuff in STD

template<class... Ts> struct tuple_size<boost::multi::detail::tuple<Ts...>> {
	static constexpr std::size_t value = sizeof...(Ts);
};

template<class T0, class... Ts> struct tuple_element<0, boost::multi::detail::tuple<T0, Ts...>> {
	using type = T0;
};

template<std::size_t N, class T0, class... Ts> struct tuple_element<N, boost::multi::detail::tuple<T0, Ts...>> {
	using type = typename tuple_element<N - 1, boost::multi::detail::tuple<Ts...>>::type;
};

template<std::size_t N, class... Ts>
constexpr auto get(boost::multi::detail::tuple<Ts...> const& t)
->decltype(boost::multi::detail::get<N>(t)) {
	return boost::multi::detail::get<N>(t); }

template<std::size_t N, class... Ts>
constexpr auto get(boost::multi::detail::tuple<Ts...> & t)
->decltype(boost::multi::detail::get<N>(t)) {
	return boost::multi::detail::get<N>(t); }

template<std::size_t N, class... Ts>
constexpr auto get(boost::multi::detail::tuple<Ts...> && t)
->decltype(boost::multi::detail::get<N>(std::move(t))) {
	return boost::multi::detail::get<N>(std::move(t)); }

template <class F, class Tuple, std::size_t... I>
constexpr auto apply_timpl(F&& f, Tuple&& t, std::index_sequence<I...> /*012*/) -> decltype(auto) {
	(void)t;  // fix "error #827: parameter "t" was never referenced" in NVC++ and "error #869: parameter "t" was never referenced" in oneAPI-ICPC
    return std::forward<F>(f)(boost::multi::detail::get<I>(std::forward<Tuple>(t))...);
}

template <class F, class... Ts>
constexpr auto apply(F&& f, boost::multi::detail::tuple<Ts...> const& t) -> decltype(auto) {
    return apply_timpl(
        std::forward<F>(f), t,
        std::make_index_sequence<sizeof...(Ts)>{}
	);
}

template <class F, class... Ts>
constexpr auto apply(F&& f, boost::multi::detail::tuple<Ts...>& t) -> decltype(auto) {
    return apply_timpl(
        std::forward<F>(f), t,
        std::make_index_sequence<sizeof...(Ts)>{}
	);
}

template <class F, class... Ts>
constexpr auto apply(F&& f, boost::multi::detail::tuple<Ts...>&& t) -> decltype(auto) {
    return apply_timpl(
        std::forward<F>(f), std::move(t),
        std::make_index_sequence<sizeof...(Ts)>{}
	);
}

}  // end namespace std

namespace boost::multi {  // NOLINT(modernize-concat-nested-namespaces) keep c++14 compat
namespace detail {
//template<class... Ts>
//constexpr auto make_tuple(Ts... ts) -> tuple<std::decay_t<Ts>...> {
//	return tuple<std::decay_t<Ts>...>(std::move(ts)...) }

//template<class Tuple1, std::size_t... Indices>
//auto tuple_zip_impl(Tuple1&& t1, std::index_sequence<Indices...> /*012*/) {
//	return make_tuple(
//		make_tuple(
//			std::get<Indices>(std::forward<Tuple1>(t1))
//		)...
//	);
//}

//template<class Tuple1, class Tuple2, std::size_t... Is>
//auto tuple_zip_impl(Tuple1&& t1, Tuple2&& t2, std::index_sequence<Is...> /*012*/) {
//	return make_tuple(
//		make_tuple(
//			std::get<Is>(std::forward<Tuple1>(t1)),
//			std::get<Is>(std::forward<Tuple2>(t2))
//		)...
//	);
//}

template<class Tuple1, class Tuple2, std::size_t... Is>
constexpr auto tuple_zip_impl(Tuple1&& t1, Tuple2&& t2, std::index_sequence<Is...> /*012*/) {
	using boost::multi::detail::get;
	return boost::multi::detail::mk_tuple(
		boost::multi::detail::mk_tuple(
			get<Is>(std::forward<Tuple1>(t1)),
			get<Is>(std::forward<Tuple2>(t2))
		)...
	);
}

template<class Tuple1, class Tuple2, class Tuple3, std::size_t... Is>
constexpr auto tuple_zip_impl(Tuple1&& t1, Tuple2&& t2, Tuple3&& t3, std::index_sequence<Is...> /*012*/) {
	using boost::multi::detail::get;
	return boost::multi::detail::mk_tuple(
		boost::multi::detail::mk_tuple(
			get<Is>(std::forward<Tuple1>(t1)),
			get<Is>(std::forward<Tuple2>(t2)),
			get<Is>(std::forward<Tuple3>(t3))
		)...
	);
}

template<class Tuple1, class Tuple2, class Tuple3, class Tuple4, std::size_t... Is>
constexpr auto tuple_zip_impl(Tuple1&& t1, Tuple2&& t2, Tuple3&& t3, Tuple4&& t4, std::index_sequence<Is...> /*012*/) {
	using boost::multi::detail::get;
	return boost::multi::detail::mk_tuple(
		boost::multi::detail::mk_tuple(
			get<Is>(std::forward<Tuple1>(t1)),
			get<Is>(std::forward<Tuple2>(t2)),
			get<Is>(std::forward<Tuple3>(t3)),
			get<Is>(std::forward<Tuple4>(t4))
		)
		...
	);
}

//template<class Tuple1, class... Tuples>
//auto tuple_zip(Tuple1&& t1, Tuples&&... ts) {
//	return detail::tuple_zip_impl(
//		std::forward<Tuple1>(t1), std::forward<Tuples>(ts)...,
//		std::make_index_sequence<std::tuple_size<typename std::decay<Tuple1>::type>::value>()
//	);
//}

template<class T1, class T2>
constexpr auto tuple_zip(T1&& t1, T2&& t2) {
	return detail::tuple_zip_impl(
		std::forward<T1>(t1), std::forward<T2>(t2),
		std::make_index_sequence<std::tuple_size<typename std::decay<T1>::type>::value>()
	);
}

template<class T1, class T2, class T3>
constexpr auto tuple_zip(T1&& t1, T2&& t2, T3&& t3) {
	return detail::tuple_zip_impl(
		std::forward<T1>(t1), std::forward<T2>(t2), std::forward<T3>(t3),
		std::make_index_sequence<std::tuple_size<typename std::decay<T1>::type>::value>()
	);
}

template<class T1, class T2, class T3, class T4>
constexpr auto tuple_zip(T1&& t1, T2&& t2, T3&& t3, T4&& t4) {
	return detail::tuple_zip_impl(
		std::forward<T1>(t1), std::forward<T2>(t2), std::forward<T3>(t3), std::forward<T4>(t4),
		std::make_index_sequence<std::tuple_size<typename std::decay<T1>::type>::value>()
	);
}


}  // end namespace detail
}  // end namespace boost::multi

#endif
