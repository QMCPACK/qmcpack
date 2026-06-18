// Copyright 2021-2026 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_MULTI_DETAIL_TUPLE_ZIP_HPP
#define BOOST_MULTI_DETAIL_TUPLE_ZIP_HPP
// #pragma once

// #include <cstddef>      // for size_t
#include <tuple>        // for deprecated functions  // for make_index_sequence, index_sequence, tuple_element, tuple_size, apply, tuple
#include <type_traits>  // for declval, decay_t, conditional_t, true_type
#include <utility>      // for forward, move

#ifdef __NVCC__
#define BOOST_MULTI_HD __host__ __device__
#else
#define BOOST_MULTI_HD
#endif

#ifdef __NVCC__
#define BOOST_MULTI_DEV __device__
#else
#define BOOST_MULTI_DEV
#endif

// namespace boost::multi 
// {}
// #ifdef __NVCC__
// template<class F>
// struct device {
// 	F f_;
// 	template<class... Ts> __device__ auto operator()(Ts&& ... ts) const {
// 		return f_(std::forward<Ts>(ts)...);
// 	}
// }
// #else
// #define BOOST_MULTI_DEV
// #endif

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable : 4514)  // boost::multi::detail::tuple<>::operator <': unreferenced inline function has been removed
#pragma warning(disable : 4623)  // default constructor was implicitly defined as deleted
#pragma warning(disable : 4626)  // assignment operator was implicitly defined as deleted
#endif

namespace boost::multi {  // NOLINT(modernize-concat-nested-namespaces) keep c++14 compat
namespace detail {

// we need a custom tuple type, so some fundamental types (e.g. iterators) are trivially constructible
template<class... Ts> class tuple;  // TODO(correaa) consider renaming it to `tpl`

template<> class tuple<> {  // NOLINT(cppcoreguidelines-special-member-functions,hicpp-special-member-functions)
 public:
	tuple()             = default;
	tuple(tuple const&) = default;

	auto operator=(tuple const&) -> tuple& = default;

	BOOST_MULTI_HD constexpr auto operator==(tuple const& /*other*/) const { return true; }
	BOOST_MULTI_HD constexpr auto operator!=(tuple const& /*other*/) const  { return false; }

	BOOST_MULTI_HD constexpr auto operator<(tuple const& /*other*/) const { return false; }
	BOOST_MULTI_HD constexpr auto operator>(tuple const& /*other*/) const { return false; }

	template<class F>
	BOOST_MULTI_HD constexpr friend auto apply(F&& fn, tuple<> const& /*self*/) -> decltype(auto) {  // NOLINT(cert-dcl58-cpp) normal idiom to defined tuple get
		return std::forward<F>(fn)();
	}
};

template<class T0, class... Ts> class tuple<T0, Ts...> : tuple<Ts...> {  // NOLINT(cppcoreguidelines-special-member-functions,hicpp-special-member-functions)
	T0 head_;  // NOLINT(cppcoreguidelines-avoid-const-or-ref-data-members) can be a reference
	using head_type = T0;
	using tail_type = tuple<Ts...>;

 public:
	BOOST_MULTI_HD constexpr auto head() const& -> T0 const& { return head_; }
	BOOST_MULTI_HD constexpr auto head() && -> T0&& { return std::move(head_); }  //-V::659 this function moves the elements, not the whole object
	BOOST_MULTI_HD constexpr auto head() & -> T0& { return head_; }

	BOOST_MULTI_HD constexpr auto tail() const& -> tail_type const& { return static_cast<tail_type const&>(*this); }
	BOOST_MULTI_HD constexpr auto tail() && -> decltype(auto) { return static_cast<tail_type&&>(*this); }
	BOOST_MULTI_HD constexpr auto tail() & -> tail_type& { return static_cast<tail_type&>(*this); }

	constexpr tuple()             = default;
	constexpr tuple(tuple const&) = default;

	// this is horrible hack and can produce ODR reported by Circle
	// operator std::tuple<T0&, Ts&...>() & {  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions,cppcoreguidelines-explicit-constructor,misc-explicit-constructor)
	// 	using std::apply;
	// 	return apply([](auto&&... ts) {return std::tuple<T0&, Ts&...>{std::forward<decltype(ts)>(ts)...}; }, *this);
	// }
	// operator std::tuple<T0&, Ts&...>() && {  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions,cppcoreguidelines-explicit-constructor,misc-explicit-constructor)
	// 	using std::apply;
	// 	return apply([](auto&&... ts) {return std::tuple<T0&, Ts&...>{std::forward<decltype(ts)>(ts)...}; }, *this);
	// }
	// operator std::tuple<T0&, Ts&...>() const& {  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions,cppcoreguidelines-explicit-constructor,misc-explicit-constructor)
	// 	using std::apply;
	// 	return apply([](auto&&... ts) {return std::tuple<T0&, Ts&...>{std::forward<decltype(ts)>(ts)...}; }, *this);
	// }

	// TODO(correaa) make conditional explicit constructor depending on the conversions for T0, Ts...
	BOOST_MULTI_HD constexpr explicit tuple(T0 head, tuple<Ts...> tail) : tail_type{std::move(tail)}, head_{std::move(head)} {}
	// cppcheck-suppress noExplicitConstructor ; allow bracket init in function argument // NOLINTNEXTLINE(runtime/explicit)
	BOOST_MULTI_HD constexpr tuple(T0 head, Ts... tail) : tail_type{tail...}, head_{head} {}  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions,cppcoreguidelines-explicit-constructor,misc-explicit-constructor) to allow bracket function calls

	// cppcheck-suppress noExplicitConstructor ; allow bracket init in function argument // NOLINTNEXTLINE(runtime/explicit)
	template<class TT0 = T0,
		std::enable_if_t<!std::is_same_v<TT0, tuple>, int> =0,  // NOLINT(modernize-use-constraints,modernize-type-traits) for C++20
		std::enable_if_t<sizeof(TT0*) && (sizeof...(Ts) == 0), int> =0  // NOLINT(modernize-use-constraints) for C++20
	>
	// cppcheck-suppress noExplicitConstructor ; see below
	BOOST_MULTI_HD constexpr tuple(TT0 head) : tail_type{}, head_{head} {}  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions,cppcoreguidelines-explicit-constructor,misc-explicit-constructor) to allow bracket function calls

	// cppcheck-suppress noExplicitConstructor ; allow bracket init in function argument // NOLINTNEXTLINE(runtime/explicit)
	BOOST_MULTI_HD constexpr explicit tuple(::std::tuple<T0, Ts...> other) : tuple(::std::apply([](auto... es) -> auto {return tuple(es...);}, other)) {}  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions,cppcoreguidelines-explicit-constructor,misc-explicit-constructor)

	constexpr auto operator=(tuple const&) -> tuple& = default;

	template<class... TTs>
	constexpr auto operator=(tuple<TTs...> const& other)  // NOLINT(cppcoreguidelines-c-copy-assignment-signature,misc-unconventional-assign-operator) signature used for SFINAE
		-> decltype(std::declval<head_type&>() = other.head(), std::declval<tail_type&>() = other.tail(), std::declval<tuple&>()) {
		head_ = other.head(), tail() = other.tail();
		return *this;
	}

	BOOST_MULTI_HD constexpr auto operator==(tuple const& other) const -> bool {
		return
			head_ == other.head_
			&& 
			tail() == other.tail()
		;
	}
	BOOST_MULTI_HD constexpr auto operator!=(tuple const& other) const -> bool { return head_ != other.head_ || tail() != other.tail(); }

	BOOST_MULTI_HD constexpr auto operator<(tuple const& other) const {
		if(head_ < other.head_) {
			return true;
		}
		if(other.head_ < head_) {
			return false;
		}
		return tail() < other.tail();
	}
	BOOST_MULTI_HD constexpr auto operator>(tuple const& other) const {
		if(head_ > other.head_) {
			return true;
		}
		if(other.head_ > head_) {
			return false;
		}
		return tail() > other.tail();
	}

 private:

#ifdef __NVCC__
#pragma nv_exec_check_disable
#endif
	template<class F, std::size_t... I>
	BOOST_MULTI_HD constexpr auto apply_impl_(F&& fn, std::index_sequence<I...> /*012*/) const& -> decltype(auto) {  // NOLINT(cert-dcl58-cpp) normal idiom to defined tuple get
		return std::forward<F>(fn)(this->get<I>()...);
	}

#ifdef __NVCC__
#pragma nv_exec_check_disable
#endif
	template<class F, std::size_t... I>
	BOOST_MULTI_DEV constexpr auto apply_impl_(F&& fn, std::index_sequence<I...> /*012*/) & -> decltype(auto) {  // NOLINT(cert-dcl58-cpp) normal idiom to defined tuple get
		return std::forward<F>(fn)(this->get<I>()...);
	}

#ifdef __NVCC__
#pragma nv_exec_check_disable
#endif
	template<class F, std::size_t... I>
	BOOST_MULTI_HD  // BOOST_MULTI_DEV
	constexpr auto apply_impl_(F&& fn, std::index_sequence<I...> /*012*/) && -> decltype(auto) {  // NOLINT(cert-dcl58-cpp) normal idiom to defined tuple get
		return std::forward<F>(fn)(std::move(*this).template get<I>()...);
	}

 public:
	template<class F>
	BOOST_MULTI_HD constexpr auto apply(F&& fn) const& -> decltype(auto) {  // NOLINT(cert-dcl58-cpp) normal idiom to defined tuple get
		return apply_impl_(std::forward<F>(fn), std::make_index_sequence<sizeof...(Ts) + 1>{});
	}
	template<class F>
	BOOST_MULTI_HD constexpr auto apply(F&& fn) & -> decltype(auto) {  // NOLINT(cert-dcl58-cpp) normal idiom to defined tuple get
		return apply_impl_(std::forward<F>(fn), std::make_index_sequence<sizeof...(Ts) + 1>{});
	}
	template<class F>
	BOOST_MULTI_HD constexpr auto apply(F&& fn) && -> decltype(auto) {  // NOLINT(cert-dcl58-cpp) normal idiom to defined tuple get
		return std::move(*this).apply_impl_(std::forward<F>(fn), std::make_index_sequence<sizeof...(Ts) + 1>{});
	}

	template<class F>
	friend BOOST_MULTI_HD constexpr auto apply(F&& fn, tuple<T0, Ts...> const& self) -> decltype(auto) {  // NOLINT(cert-dcl58-cpp) normal idiom to defined tuple get
		return self.apply(std::forward<F>(fn));
	}

	template<class F>
	friend BOOST_MULTI_HD constexpr auto apply(F&& fn, tuple<T0, Ts...> & self) -> decltype(auto) {  // NOLINT(cert-dcl58-cpp) normal idiom to defined tuple get
		return self.apply(std::forward<F>(fn));
	}

	template<class F>
	friend BOOST_MULTI_HD constexpr auto apply(F&& fn, tuple<T0, Ts...> && self) -> decltype(auto) {  // NOLINT(cert-dcl58-cpp) normal idiom to defined tuple get
		return std::move(self).apply(std::forward<F>(fn));
	}

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

	template<std::size_t N, std::enable_if_t<(N==0), int> =0>  // NOLINT(modernize-use-constraints) for C++20
	BOOST_MULTI_HD constexpr auto get() const& -> T0 const& {  // NOLINT(readability-identifier-length) std naming
		return head();
	}

	template<std::size_t N, std::enable_if_t<(N!=0), int> =0>  // NOLINT(modernize-use-constraints) for C++20
	BOOST_MULTI_HD constexpr auto get() const& -> auto const& {  // NOLINT(readability-identifier-length) std naming
		return this->tail().template get<N - 1>();  // this-> for msvc 19.14 compilation
	}


#ifdef __NVCC__
	#pragma nv_diagnostic push
	#pragma nv_diag_suppress = implicit_return_from_non_void_function  // in place of global -Xcudafe \"--diag_suppress=implicit_return_from_non_void_function\"
#elif defined(__NVCOMPILER)  // NOLINT(readability-use-concise-preprocessor-directives) for C++23
	#pragma diagnostic push
	#pragma diag_suppress = implicit_return_from_non_void_function
#elif (defined(__GNUC__) && !defined(__EDG__))
	#pragma GCC diagnostic push
	#pragma GCC diagnostic ignored "-Wreturn-type"
#endif

	template<std::size_t N>
	BOOST_MULTI_HD constexpr auto get() & -> decltype(auto) {  // NOLINT(readability-identifier-length) std naming
		if constexpr(N == 0) {
			return head();
		} else {
			return tail().template get<N - 1>();
		}
	}

	template<std::size_t N>
	BOOST_MULTI_HD constexpr auto get() && -> decltype(auto) {  // NOLINT(readability-identifier-length) std naming
		if constexpr(N == 0) {
			return std::move(*this).head();
		} else {
			return std::move(*this).tail().template get<N - 1>();
		}
	}
};

#ifdef __NVCC__
#pragma nv_diagnostic pop
#elif defined(__NVCOMPILER)  // NOLINT(readability-use-concise-preprocessor-directives) for C++23
#pragma diagnostic pop
#elif (defined(__GNUC__) && !defined(__EDG__))
#pragma GCC diagnostic pop
#endif

#ifdef _MSC_VER
#pragma warning(pop)
#endif

#ifdef __INTEL_COMPILER  // this instance is necessary due to a bug in intel compiler icpc
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

#if defined(__cpp_deduction_guides) && (__cpp_deduction_guides >= 201703L) 
template<class T0, class... Ts> tuple(T0, tuple<Ts...>) -> tuple<T0, Ts...>;
#endif

template<class T0, class... Ts> constexpr auto mk_tuple(T0 head, Ts... tail) {
	return tuple<T0, Ts...>(std::move(head), std::move(tail)...);
}

template<class T0, class... Ts> constexpr auto tie(T0& head, Ts&... tail) {
	return tuple<T0&, Ts&...>(head, tail...);
}

template<class T0, class... Ts> BOOST_MULTI_HD constexpr auto ht_tuple(T0 head, tuple<Ts...> tail)
-> tuple<T0, Ts...> {
	return tuple<T0, Ts...>(std::move(head), std::move(tail));
}

template<class T0, class Tuple> struct tuple_prepend;

template<class T0, class... Ts>
struct tuple_prepend<T0, tuple<Ts...>> {
	using type = tuple<T0, Ts...>;
};

template<class T0, class Tuple>
using tuple_prepend_t = typename tuple_prepend<T0, Tuple>::type;  // NOLINT(readability-redundant-typename) for C++20

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
BOOST_MULTI_HD constexpr auto tail(tuple<T0, Ts...> const& t) -> decltype(t.tail()) { return t.tail(); }  // NOLINT(readability-identifier-length) std naming

template<class T0, class... Ts>
BOOST_MULTI_HD constexpr auto tail(tuple<T0, Ts...>&& t) -> decltype(std::move(t).tail()) { return std::move(t).tail(); }  // NOLINT(readability-identifier-length) std naming

template<class T0, class... Ts>
BOOST_MULTI_HD constexpr auto tail(tuple<T0, Ts...>& t) -> decltype(t.tail()) { return t.tail(); }  // NOLINT(readability-identifier-length) std naming

#ifdef __NVCC__
	#pragma nv_diagnostic push
	#pragma nv_diag_suppress = implicit_return_from_non_void_function  // in place of global -Xcudafe \"--diag_suppress=implicit_return_from_non_void_function\"
#elif defined(__NVCOMPILER)  // NOLINT(readability-use-concise-preprocessor-directives) for C++23
	#pragma diagnostic push
	#pragma diag_suppress = implicit_return_from_non_void_function
#elif defined(__GNUC__) && !defined(__EDG__)
	#pragma GCC diagnostic push
	#pragma GCC diagnostic ignored "-Wreturn-type"
#endif

template<std::size_t N, class T0, class... Ts>
BOOST_MULTI_HD constexpr auto get(tuple<T0, Ts...> const& t) -> auto const& {  // NOLINT(readability-identifier-length) std naming
	if constexpr(N == 0) {
		return t.head();
	} else {
		using std::get;
		return get<N - 1>(t.tail());
	}
}

template<std::size_t N, class T0, class... Ts>
BOOST_MULTI_HD constexpr auto get(tuple<T0, Ts...>& tup) -> auto& {
	if constexpr(N == 0) {
		return tup.head();
	} else {
		return get<N - 1>(tup.tail());
	}
}

template<std::size_t N, class T0, class... Ts>
BOOST_MULTI_HD constexpr auto get(tuple<T0, Ts...>&& tup) -> auto&& {  //-V::659 overload for move
	if constexpr(N == 0) {
		return std::move(std::move(tup)).head();
	} else {
		return get<N - 1>(std::move(std::move(tup).tail()));
	}
}

#ifdef __NVCC__
	#pragma nv_diagnostic pop
#elif defined(__NVCOMPILER)  // NOLINT(readability-use-concise-preprocessor-directives) for C++23
	#pragma diagnostic pop
#elif defined(__GNUC__) && !defined(__EDG__)
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
struct std::tuple_size<boost::multi::detail::tuple<Ts...>> : std::integral_constant<std::size_t, sizeof...(Ts)> {  // NOLINT(cert-dcl58-cpp,bugprone-std-namespace-modification) for structured bindings
	// // cppcheck-suppress unusedStructMember
	// static constexpr std::size_t value = sizeof...(Ts);
};

template<>
struct std::tuple_element<0, boost::multi::detail::tuple<>> {  // NOLINT(cert-dcl58-cpp) to have structured bindings
	using type = void;
};

template<class T0>
struct std::tuple_element<0, boost::multi::detail::tuple<T0>> {  // NOLINT(cert-dcl58-cpp,bugprone-std-namespace-modification) for structured bindings
	using type = T0;
};

template<class T0, class... Ts>
struct std::tuple_element<0, boost::multi::detail::tuple<T0, Ts...>> {  // NOLINT(cert-dcl58-cpp,bugprone-std-namespace-modification) for structured bindings
	using type = T0;
};

template<std::size_t N, class T0, class... Ts>
struct std::tuple_element<N, boost::multi::detail::tuple<T0, Ts...>> {  // NOLINT(cert-dcl58-cpp,bugprone-std-namespace-modification) for structured bindings
	using type = tuple_element_t<N - 1, boost::multi::detail::tuple<Ts...>>;
};

#ifdef __NVCC__
#pragma nv_exec_check_disable
#endif
template<class F, class Tuple, std::size_t... I>
BOOST_MULTI_HD constexpr auto std_apply_timpl(F&& fn, Tuple&& tp, std::index_sequence<I...> /*012*/) -> decltype(auto) {  // NOLINT(cert-dcl58-cpp) normal idiom to defined tuple get
	(void)tp;  // fix "error #827: parameter "t" was never referenced" in NVC++ and "error #869: parameter "t" was never referenced" in oneAPI-ICPC
	return std::forward<F>(fn)(boost::multi::detail::get<I>(std::forward<Tuple>(tp))...);  // NOLINT(bugprone-use-after-move,hicpp-invalid-access-moved) use forward_as?
}

namespace std {  // NOLINT(cert-dcl58-cpp,bugprone-std-namespace-modification) to implement structured bindings

template<class F, class... Ts>
BOOST_MULTI_HD constexpr auto apply(F&& fn, boost::multi::detail::tuple<Ts...> const& tp) -> decltype(auto) {  // NOLINT(cert-dcl58-cpp,bugprone-std-namespace-modification) normal to define tuple get
	return std_apply_timpl(
		std::forward<F>(fn), tp,
		std::make_index_sequence<sizeof...(Ts)>{}
	);
}

template<class F, class... Ts>
BOOST_MULTI_HD constexpr auto apply(F&& fn, boost::multi::detail::tuple<Ts...>& tp) -> decltype(auto) {  // NOLINT(cert-dcl58-cpp,bugprone-std-namespace-modification) normal to define tuple get
	return std_apply_timpl(
		std::forward<F>(fn), tp,
		std::make_index_sequence<sizeof...(Ts)>{}
	);
}

template<class F, class... Ts>
BOOST_MULTI_HD constexpr auto apply(F&& fn, boost::multi::detail::tuple<Ts...>&& tp) -> decltype(auto) {  // NOLINT(cert-dcl58-cpp,bugprone-std-namespace-modification) normal idiom to defined tuple get
	return std_apply_timpl(	
		std::forward<F>(fn), std::move(tp),
		std::make_index_sequence<sizeof...(Ts)>{}
	);
}

}  // end namespace std

#ifdef __clang__
#pragma clang diagnostic pop
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
			get<Is>(std::forward<Tuple1>(tup1)),  // NOLINT(bugprone-use-after-move,hicpp-invalid-access-moved) use forward_as
			get<Is>(std::forward<Tuple2>(tup2)),  // NOLINT(bugprone-use-after-move,hicpp-invalid-access-moved) use forward_as
			get<Is>(std::forward<Tuple3>(tup3)),  // NOLINT(bugprone-use-after-move,hicpp-invalid-access-moved) use forward_as
			get<Is>(std::forward<Tuple4>(tup4))   // NOLINT(bugprone-use-after-move,hicpp-invalid-access-moved) use forward_as
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
		std::make_index_sequence<std::tuple_size_v<std::decay_t<T1>>>()
	);
}

template<class T1, class T2, class T3>
constexpr auto tuple_zip(T1&& tup1, T2&& tup2, T3&& tup3) {
	return detail::tuple_zip_impl(
		std::forward<T1>(tup1), std::forward<T2>(tup2), std::forward<T3>(tup3),
		std::make_index_sequence<std::tuple_size_v<std::decay_t<T1>>>()
	);
}

template<class T1, class T2, class T3, class T4>
constexpr auto tuple_zip(T1&& tup1, T2&& tup2, T3&& tup3, T4&& tup4) {
	return detail::tuple_zip_impl(
		std::forward<T1>(tup1), std::forward<T2>(tup2), std::forward<T3>(tup3), std::forward<T4>(tup4),
		std::make_index_sequence<std::tuple_size_v<std::decay_t<T1>>>()
	);
}

template<class T1, class T2, class T3, class T4, class T5>
constexpr auto tuple_zip(T1&& tup1, T2&& tup2, T3&& tup3, T4&& tup4, T5&& tup5) {
	return detail::tuple_zip_impl(
		std::forward<T1>(tup1), std::forward<T2>(tup2), std::forward<T3>(tup3), std::forward<T4>(tup4), std::forward<T5>(tup5),
		std::make_index_sequence<std::tuple_size_v<std::decay_t<T1>>>()
	);
}

template<class T, class Tuple, class = std::make_index_sequence<std::tuple_size_v<Tuple>>>
struct all_elements_convertible_to;

template<class T, class Tuple, std::size_t... Is>
struct all_elements_convertible_to<T, Tuple, std::index_sequence<Is...>>
: std::conjunction<std::is_convertible<std::tuple_element_t<Is, Tuple>, T>...> {};

}  // end namespace detail

using detail::tie;

}  // end namespace boost::multi

#undef BOOST_MULTI_HD

#endif  // BOOST_MULTI_DETAIL_TUPLE_ZIP_HPP
