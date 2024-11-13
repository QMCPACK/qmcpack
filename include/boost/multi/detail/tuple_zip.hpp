// Copyright 2021-2024 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_MULTI_DETAIL_TUPLE_ZIP_HPP
#define BOOST_MULTI_DETAIL_TUPLE_ZIP_HPP
#pragma once

#include <cstddef>      // for size_t
#include <tuple>        // for deprecated functions  // for make_index_sequence, index_sequence, tuple_element, tuple_size, apply, tuple
#include <type_traits>  // for declval, decay_t, conditional_t, true_type
#include <utility>      // for forward, move

namespace boost::multi {  // NOLINT(modernize-concat-nested-namespaces) keep c++14 compat
namespace detail {

template<class... Ts> class tuple;

template<> class tuple<> {  // NOLINT(cppcoreguidelines-special-member-functions,hicpp-special-member-functions)
 public:
	constexpr tuple()             = default;
	constexpr tuple(tuple const&) = default;

	constexpr auto operator=(tuple const&) -> tuple& = default;

	constexpr auto operator==(tuple const& /*other*/) const -> bool { return true; }
	constexpr auto operator!=(tuple const& /*other*/) const -> bool { return false; }

	constexpr auto operator<(tuple const& /*other*/) const { return false; }
	constexpr auto operator>(tuple const& /*other*/) const { return false; }
};

template<class T0, class... Ts> class tuple<T0, Ts...> : tuple<Ts...> {  // NOLINT(cppcoreguidelines-special-member-functions,hicpp-special-member-functions)
	T0 head_;  // NOLINT(cppcoreguidelines-avoid-const-or-ref-data-members) can be a reference
	using head_type = T0;
	using tail_type = tuple<Ts...>;

 public:
	constexpr auto head() const& -> T0 const& { return head_; }
	constexpr auto head() && -> T0&& { return std::move(head_); }
	constexpr auto head() & -> T0& { return head_; }

	constexpr auto tail() const& -> tail_type const& { return static_cast<tail_type const&>(*this); }
	constexpr auto tail() && -> decltype(auto) { return static_cast<tail_type&&>(*this); }
	constexpr auto tail() & -> tail_type& { return static_cast<tail_type&>(*this); }

	constexpr tuple()             = default;
	constexpr tuple(tuple const&) = default;

	// TODO(correaa) make conditional explicit constructor depending on the conversions for T0, Ts...
	constexpr explicit tuple(T0 head, tuple<Ts...> tail) : tail_type{std::move(tail)}, head_{std::move(head)} {}
	// cppcheck-suppress noExplicitConstructor ; allow bracket init in function argument // NOLINTNEXTLINE(runtime/explicit)
	constexpr tuple(T0 head, Ts... tail) : tail_type{tail...}, head_{head} {}  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions) to allow bracket function calls

	// cppcheck-suppress noExplicitConstructor ; allow bracket init in function argument // NOLINTNEXTLINE(runtime/explicit)
	constexpr explicit tuple(::std::tuple<T0, Ts...> other) : tuple(::std::apply([](auto... es) {return tuple(es...);}, other)) {}  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions)

	constexpr auto operator=(tuple const&) -> tuple& = default;

	template<class... TTs>
	constexpr auto operator=(tuple<TTs...> const& other)  // NOLINT(cppcoreguidelines-c-copy-assignment-signature,misc-unconventional-assign-operator) signature used for SFINAE
		-> decltype(std::declval<head_type&>() = other.head(), std::declval<tail_type&>() = other.tail(), std::declval<tuple&>()) {
		head_ = other.head(), tail() = other.tail();
		return *this;
	}

	constexpr auto operator==(tuple const& other) const -> bool { return head_ == other.head_ && tail() == other.tail(); }
	constexpr auto operator!=(tuple const& other) const -> bool { return head_ != other.head_ || tail() != other.tail(); }

	constexpr auto operator<(tuple const& other) const {
		if(head_ < other.head_) {
			return true;
		}
		if(other.head_ < head_) {
			return false;
		}
		return tail() < other.tail();
	}
	constexpr auto operator>(tuple const& other) const {
		if(head_ > other.head_) {
			return true;
		}
		if(other.head_ > head_) {
			return false;
		}
		return tail() > other.tail();
	}

 private:

	template<class F, std::size_t... I>
	constexpr auto apply_impl_(F&& fn, std::index_sequence<I...> /*012*/) const -> decltype(auto) {  // NOLINT(cert-dcl58-cpp) normal idiom to defined tuple get
		return std::forward<F>(fn)(this->get<I>()...);
	}

 public:
	template<class F>
	constexpr auto apply(F&& fn) const -> decltype(auto) {  // NOLINT(cert-dcl58-cpp) normal idiom to defined tuple get
		return apply_impl_(std::forward<F>(fn), std::make_index_sequence<sizeof...(Ts) + 1>{});
	}

	// constexpr explicit operator std::tuple<T0, Ts...>() const {return this->apply([](auto const&... xs) {return std::tuple<T0, Ts...>(xs...);});}

 private:
	template<std::size_t N> struct priority : std::conditional_t<N == 0, std::true_type, priority<N - 1>> {};

	template<class Index>
	constexpr auto at_aux_(priority<0> /*prio*/, Index idx) const
		-> decltype(ht_tuple(std::declval<head_type const&>(), std::declval<tail_type const&>()[idx])) {
		return ht_tuple(head(), tail()[idx]);
	}

	template<class Index>
	constexpr auto at_aux_(priority<1> /*prio*/, Index idx) const
		-> decltype(ht_tuple(std::declval<head_type const&>()[idx], std::declval<tail_type const&>())) {
		return ht_tuple(head()[idx], tail());
	}

 public:
	template<class Index>
	constexpr auto operator[](Index idx) const
		-> decltype(std::declval<tuple<T0, Ts...> const&>().at_aux_(priority<1>{}, idx)) {
		return this->at_aux_(priority<1>{}, idx);
	}

	template<std::size_t N>
	constexpr auto get() const& -> auto const& {  // NOLINT(readability-identifier-length) std naming
		if constexpr(N == 0) {
			return head();
		} else {
			return tail().template get<N - 1>();
		}
	}

	template<std::size_t N>
	constexpr auto get() & -> decltype(auto) {  // NOLINT(readability-identifier-length) std naming
		if constexpr(N == 0) {
			return head();
		} else {
			return tail().template get<N - 1>();
		}
	}

	template<std::size_t N>
	constexpr auto get() && -> decltype(auto) {  // NOLINT(readability-identifier-length) std naming
		if constexpr(N == 0) {
			return std::move(*this).head();
		} else {
			return std::move(*this).tail().template get<N - 1>();
		}
	}
};

#if defined(__INTEL_COMPILER)  // this instance is necessary due to a bug in intel compiler icpc
//  TODO(correaa) : this class can be collapsed with the general case with [[no_unique_address]] in C++20
template<class T0> class tuple<T0> {  // NOLINT(cppcoreguidelines-special-member-functions,hicpp-special-member-functions)
	T0      head_;
	tuple<> tail_;

 public:
	constexpr auto head() const& -> T0 const& { return head_; }
	constexpr auto head() && -> T0&& { return std::move(head_); }
	constexpr auto head() & -> T0& { return head_; }

	constexpr auto tail() const& -> tuple<> const& { return tail_; }
	constexpr auto tail() && -> tuple<>&& { return std::move(tail_); }
	constexpr auto tail() & -> tuple<>& { return tail_; }

	constexpr tuple()             = default;
	constexpr tuple(tuple const&) = default;  // cppcheck-suppress noExplicitConstructor ; workaround cppcheck 2.11

	// cppcheck-suppress noExplicitConstructor ; allow bracket init in function argument // NOLINTNEXTLINE(runtime/explicit)
	constexpr tuple(T0 t0, tuple<> sub) : head_{std::move(t0)}, tail_{sub} {}
	constexpr explicit tuple(T0 t0) : head_{std::move(t0)}, tail_{} {}

	constexpr auto operator=(tuple const& other) -> tuple& = default;

	constexpr auto operator==(tuple const& other) const { return head_ == other.head_; }
	constexpr auto operator!=(tuple const& other) const { return head_ != other.head_; }

	constexpr auto operator<(tuple const& other) const { return head_ < other.head_; }
	constexpr auto operator>(tuple const& other) const { return head_ > other.head_; }
};
#endif

template<class T0, class... Ts> tuple(T0, tuple<Ts...>) -> tuple<T0, Ts...>;

template<class T0, class... Ts> constexpr auto mk_tuple(T0 head, Ts... tail) {
	return tuple<T0, Ts...>(std::move(head), std::move(tail)...);
}

template<class T0, class... Ts> constexpr auto tie(T0& head, Ts&... tail) {
	return tuple<T0&, Ts&...>(head, tail...);
}

template<class T0, class... Ts> constexpr auto ht_tuple(T0 head, tuple<Ts...> tail) {
	return tuple<T0, Ts...>(std::move(head), std::move(tail));
}

template<class T0, class Tuple> struct tuple_prepend;

template<class T0, class... Ts>
struct tuple_prepend<T0, tuple<Ts...>> {
	using type = tuple<T0, Ts...>;
};

template<class T0, class Tuple>
using tuple_prepend_t = typename tuple_prepend<T0, Tuple>::type;

template<class T0, class... Ts>
constexpr auto head(tuple<T0, Ts...> const& t) -> decltype(auto) {  // NOLINT(readability-identifier-length) std naming
	return t.head();
}

template<class T0, class... Ts>
constexpr auto head(tuple<T0, Ts...>&& t) -> decltype(auto) {  // NOLINT(readability-identifier-length) std naming
	return std::move(t).head();
}

template<class T0, class... Ts>
constexpr auto head(tuple<T0, Ts...>& t) -> decltype(auto) {  // NOLINT(readability-identifier-length) std naming
	return t.head();
}

template<class T0, class... Ts>
constexpr auto tail(tuple<T0, Ts...> const& t) -> decltype(t.tail()) { return t.tail(); }  // NOLINT(readability-identifier-length) std naming

template<class T0, class... Ts>
constexpr auto tail(tuple<T0, Ts...>&& t) -> decltype(std::move(t).tail()) { return std::move(t).tail(); }  // NOLINT(readability-identifier-length) std naming

template<class T0, class... Ts>
constexpr auto tail(tuple<T0, Ts...>& t) -> decltype(t.tail()) { return t.tail(); }  // NOLINT(readability-identifier-length) std naming

#if defined __NVCC__  // in place of global -Xcudafe \"--diag_suppress=implicit_return_from_non_void_function\"
#ifdef __NVCC_DIAG_PRAGMA_SUPPORT__
#pragma nv_diagnostic push
#pragma nv_diag_suppress = implicit_return_from_non_void_function
#else
#pragma diagnostic push
#pragma diag_suppress = implicit_return_from_non_void_function
#endif
#elif defined      __NVCOMPILER
#pragma diagnostic push
#pragma diag_suppress = implicit_return_from_non_void_function
#endif
#if ! defined(_MSC_VER)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wreturn-type"
#endif
template<std::size_t N, class T0, class... Ts>
constexpr auto get(tuple<T0, Ts...> const& t) -> auto const& {  // NOLINT(readability-identifier-length) std naming
	if constexpr(N == 0) {
		return t.head();
	} else {
		return get<N - 1>(t.tail());
	}
}
#if defined __NVCC__
#ifdef __NVCC_DIAG_PRAGMA_SUPPORT__
#pragma nv_diagnostic pop
#else
#pragma diagnostic pop
#endif
#elif defined      __NVCOMPILER
#pragma diagnostic pop
#endif

template<std::size_t N, class T0, class... Ts>
constexpr auto get(tuple<T0, Ts...>& tup) -> auto& {
	if constexpr(N == 0) {
		return tup.head();
	} else {
		return get<N - 1>(tup.tail());
	}
}

template<std::size_t N, class T0, class... Ts>
constexpr auto get(tuple<T0, Ts...>&& tup) -> auto&& {
	if constexpr(N == 0) {
		return std::move(std::move(tup)).head();
	} else {
		return get<N - 1>(std::move(std::move(tup).tail()));
	}
}
#if ! defined(_MSC_VER)
#pragma GCC diagnostic pop
#endif

}  // end namespace detail
}  // end namespace boost::multi

// Some versions of Clang throw warnings that stl uses class std::tuple_size instead
// of struct std::tuple_size like it should be
#ifdef __clang__
#  pragma clang diagnostic push
#  pragma clang diagnostic ignored "-Wmismatched-tags"
#endif

template<class... Ts>
struct std::tuple_size<boost::multi::detail::tuple<Ts...>> : std::integral_constant<std::size_t, sizeof...(Ts)> {  // NOLINT(cert-dcl58-cpp) to have structured bindings
	// // cppcheck-suppress unusedStructMember
	// static constexpr std::size_t value = sizeof...(Ts);
};

template<>
struct std::tuple_element<0, boost::multi::detail::tuple<>> {  // NOLINT(cert-dcl58-cpp) to have structured bindings
	using type = void;
};

template<class T0>
struct std::tuple_element<0, boost::multi::detail::tuple<T0>> {  // NOLINT(cert-dcl58-cpp) to have structured bindings
	using type = T0;
};

template<class T0, class... Ts>
struct std::tuple_element<0, boost::multi::detail::tuple<T0, Ts...>> {  // NOLINT(cert-dcl58-cpp) to have structured bindings
	using type = T0;
};

template<std::size_t N, class T0, class... Ts>
struct std::tuple_element<N, boost::multi::detail::tuple<T0, Ts...>> {  // NOLINT(cert-dcl58-cpp) to have structured bindings
	using type = typename tuple_element<N - 1, boost::multi::detail::tuple<Ts...>>::type;
};

// namespace std {  // NOLINT(cert-dcl58-cpp) to implement structured bindings

// template<std::size_t N, class... Ts>
// constexpr auto get(boost::multi::detail::tuple<Ts...> const& tp)  // NOLINT(cert-dcl58-cpp) normal idiom to defined tuple get
//  -> decltype(boost::multi::detail::get<N>(tp)) {
//  return boost::multi::detail::get<N>(tp);
// }

// template<std::size_t N, class... Ts>
// constexpr auto get(boost::multi::detail::tuple<Ts...>& tp)  // NOLINT(cert-dcl58-cpp) normal idiom to defined tuple get
//  -> decltype(boost::multi::detail::get<N>(tp)) {
//  return boost::multi::detail::get<N>(tp);
// }

// template<std::size_t N, class... Ts>
// constexpr auto get(boost::multi::detail::tuple<Ts...>&& tp)  // NOLINT(cert-dcl58-cpp) normal idiom to defined tuple get
//  -> decltype(boost::multi::detail::get<N>(std::move(tp))) {
//  return boost::multi::detail::get<N>(std::move(tp));
// }

// }  // end namespace std

namespace std {  // NOLINT(cert-dcl58-cpp) to implement structured bindings

template<class F, class Tuple, std::size_t... I>
constexpr auto std_apply_timpl(F&& fn, Tuple&& tp, std::index_sequence<I...> /*012*/) -> decltype(auto) {  // NOLINT(cert-dcl58-cpp) normal idiom to defined tuple get
	(void)tp;  // fix "error #827: parameter "t" was never referenced" in NVC++ and "error #869: parameter "t" was never referenced" in oneAPI-ICPC
	return std::forward<F>(fn)(boost::multi::detail::get<I>(std::forward<Tuple>(tp))...);
}

template<class F, class... Ts>
constexpr auto apply(F&& fn, boost::multi::detail::tuple<Ts...> const& tp) -> decltype(auto) {  // NOLINT(cert-dcl58-cpp) normal idiom to defined tuple get
	return std_apply_timpl(
		std::forward<F>(fn), tp,
		std::make_index_sequence<sizeof...(Ts)>{}
	);
}

template<class F, class... Ts>
constexpr auto apply(F&& fn, boost::multi::detail::tuple<Ts...>& tp) -> decltype(auto) {  // NOLINT(cert-dcl58-cpp) normal idiom to defined tuple get
	return std_apply_timpl(
		std::forward<F>(fn), tp,
		std::make_index_sequence<sizeof...(Ts)>{}
	);
}

template<class F, class... Ts>
constexpr auto apply(F&& fn, boost::multi::detail::tuple<Ts...>&& tp) -> decltype(auto) {  // NOLINT(cert-dcl58-cpp) normal idiom to defined tuple get
	return std_apply_timpl(
		std::forward<F>(fn), std::move(tp),
		std::make_index_sequence<sizeof...(Ts)>{}
	);
}

}  // end namespace std

#ifdef __clang__
#  pragma clang diagnostic pop
#endif

namespace boost::multi {  // NOLINT(modernize-concat-nested-namespaces) keep c++14 compat
namespace detail {

template<class Tuple1, class Tuple2, std::size_t... Is>
constexpr auto tuple_zip_impl(Tuple1&& tup1, Tuple2&& tup2, std::index_sequence<Is...> /*012*/) {
	using boost::multi::detail::get;
	return boost::multi::detail::mk_tuple(
		boost::multi::detail::mk_tuple(
			get<Is>(std::forward<Tuple1>(tup1)),
			get<Is>(std::forward<Tuple2>(tup2))
		)...
	);
}

template<class Tuple1, class Tuple2, class Tuple3, std::size_t... Is>
constexpr auto tuple_zip_impl(Tuple1&& tup1, Tuple2&& tup2, Tuple3&& tup3, std::index_sequence<Is...> /*012*/) {
	using boost::multi::detail::get;
	return boost::multi::detail::mk_tuple(
		boost::multi::detail::mk_tuple(
			get<Is>(std::forward<Tuple1>(tup1)),
			get<Is>(std::forward<Tuple2>(tup2)),
			get<Is>(std::forward<Tuple3>(tup3))
		)...
	);
}

template<class Tuple1, class Tuple2, class Tuple3, class Tuple4, std::size_t... Is>
constexpr auto tuple_zip_impl(Tuple1&& tup1, Tuple2&& tup2, Tuple3&& tup3, Tuple4&& tup4, std::index_sequence<Is...> /*012*/) {
	using boost::multi::detail::get;
	return boost::multi::detail::mk_tuple(
		boost::multi::detail::mk_tuple(
			get<Is>(std::forward<Tuple1>(tup1)),
			get<Is>(std::forward<Tuple2>(tup2)),
			get<Is>(std::forward<Tuple3>(tup3)),
			get<Is>(std::forward<Tuple4>(tup4))
		)...
	);
}

template<class Tuple1, class Tuple2, class Tuple3, class Tuple4, class Tuple5, std::size_t... Is>
constexpr auto tuple_zip_impl(Tuple1&& tup1, Tuple2&& tup2, Tuple3&& tup3, Tuple4&& tup4, Tuple5&& tup5, std::index_sequence<Is...> /*012*/) {
	using boost::multi::detail::get;
	return boost::multi::detail::mk_tuple(
		boost::multi::detail::mk_tuple(
			get<Is>(std::forward<Tuple1>(tup1)),
			get<Is>(std::forward<Tuple2>(tup2)),
			get<Is>(std::forward<Tuple3>(tup3)),
			get<Is>(std::forward<Tuple4>(tup4)),
			get<Is>(std::forward<Tuple4>(tup5))
		)...
	);
}

template<class T1, class T2>
constexpr auto tuple_zip(T1&& tup1, T2&& tup2) {
	return detail::tuple_zip_impl(
		std::forward<T1>(tup1), std::forward<T2>(tup2),
		std::make_index_sequence<std::tuple_size<std::decay_t<T1>>::value>()
	);
}

template<class T1, class T2, class T3>
constexpr auto tuple_zip(T1&& tup1, T2&& tup2, T3&& tup3) {
	return detail::tuple_zip_impl(
		std::forward<T1>(tup1), std::forward<T2>(tup2), std::forward<T3>(tup3),
		std::make_index_sequence<std::tuple_size<std::decay_t<T1>>::value>()
	);
}

template<class T1, class T2, class T3, class T4>
constexpr auto tuple_zip(T1&& tup1, T2&& tup2, T3&& tup3, T4&& tup4) {
	return detail::tuple_zip_impl(
		std::forward<T1>(tup1), std::forward<T2>(tup2), std::forward<T3>(tup3), std::forward<T4>(tup4),
		std::make_index_sequence<std::tuple_size<std::decay_t<T1>>::value>()
	);
}

template<class T1, class T2, class T3, class T4, class T5>
constexpr auto tuple_zip(T1&& tup1, T2&& tup2, T3&& tup3, T4&& tup4, T5&& tup5) {
	return detail::tuple_zip_impl(
		std::forward<T1>(tup1), std::forward<T2>(tup2), std::forward<T3>(tup3), std::forward<T4>(tup4), std::forward<T5>(tup5),
		std::make_index_sequence<std::tuple_size<std::decay_t<T1>>::value>()
	);
}

}  // end namespace detail

using detail::tie;

}  // end namespace boost::multi
#endif
