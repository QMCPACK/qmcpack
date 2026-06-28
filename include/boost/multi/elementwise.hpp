// Copyright 2025-2026 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_MULTI_ELEMENTWISE_HPP
#define BOOST_MULTI_ELEMENTWISE_HPP

#include "boost/multi/array_ref.hpp"
#include "boost/multi/restriction.hpp"
#include "boost/multi/utility.hpp"  // for multi::detail::apply_square

#include <type_traits>

#ifdef __NVCC__
#define BOOST_MULTI_HD __host__ __device__
#else
#define BOOST_MULTI_HD
#endif

namespace boost::multi {

template<class A> struct bind_category {
	using type = A;
};

template<class T, dimensionality_type D, class... Ts>
struct bind_category<::boost::multi::array<T, D, Ts...>> {
	using type = ::boost::multi::array<T, D, Ts...>;
};

template<class T, dimensionality_type D, class... Ts>
struct bind_category<::boost::multi::array<T, D, Ts...>&> {
	using type = ::boost::multi::array<T, D, Ts...>&;
};

template<class T, dimensionality_type D, class... Ts>
struct bind_category<::boost::multi::array<T, D, Ts...> const&> {
	using type = ::boost::multi::array<T, D, Ts...> const&;
};

template<class T, dimensionality_type D, class... Ts>
struct bind_category<::boost::multi::subarray<T, D, Ts...>> {
	using type = ::boost::multi::subarray<T, D, Ts...>;
};

template<class T, dimensionality_type D, class... Ts>
struct bind_category<::boost::multi::subarray<T, D, Ts...>&> {
	using type = ::boost::multi::subarray<T, D, Ts...>&;
};

template<class T, dimensionality_type D, class... Ts>
struct bind_category<::boost::multi::subarray<T, D, Ts...> const&> {
	using type = ::boost::multi::subarray<T, D, Ts...> const&;
};

namespace elementwise {

template<class F, class A, class... Arrays, typename = decltype(std::declval<F&&>()(std::declval<typename std::decay_t<A>::reference>(), std::declval<typename std::decay_t<Arrays>::reference>()...))>
constexpr auto apply_front(F&& fun, A&& arr, Arrays&&... arrs) {
	return [fun_ = std::forward<F>(fun), &arr, &arrs...](auto is) -> decltype(auto) { return fun_(arr[is], arrs[is]...); } ^ multi::extents_t<1>({arr.extent()});
}

template<class F, class... A> struct apply_bind_t;

template<class F, class A>
struct apply_bind_t<F, A> {
	F fun_;
	A a_;

	template<class... Is>
	constexpr auto operator()(Is... is) const {
		return fun_(detail::invoke_square(a_, is...));  // a_[is...] in. C++23
	}
};

template<class F, class A, class B>
struct apply_bind_t<F, A, B> {
	F fun_;
	A a_;
	B b_;

	template<class... Is>
	constexpr auto operator()(Is... is) const {
		return fun_(
			multi::detail::invoke_square(a_, is...),  // a_[is...] in C++23
			multi::detail::invoke_square(b_, is...)   // b_[is...] in C++23
		);
	}
};

template<class F, class A, class... As, typename = decltype(std::declval<F&&>()(std::declval<typename std::decay_t<A>::element>(), std::declval<typename std::decay_t<As>::element>()...))>
constexpr auto apply(F&& fun, A&& arr, As&&... arrs) {
	auto const xs = arr.extents();  // TODO(correaa) consider storing home() cursor only
	assert(((xs == arrs.extents()) && ...));
	return multi::restricted(apply_bind_t<F, std::decay_t<A>, std::decay_t<As>...>{std::forward<F>(fun), std::forward<A>(arr), std::forward<As>(arrs)...}, xs);
}

template<class T>
class identity_bind {
	T val_;

 public:
	template<class TT, std::enable_if_t<!std::is_base_of_v<identity_bind, std::decay_t<TT>>, int> = 0>  // NOLINT(modernize-use-constraints) for C++20
	explicit constexpr identity_bind(TT&& val) : val_{std::forward<TT>(val)} {}                         // NOLINT(bugprone-forwarding-reference-overload)

	BOOST_MULTI_HD constexpr auto operator()() const -> auto& { return val_; }
};

template<class F, class A, class B>
constexpr auto map(F&& fun, A&& alpha, B&& omega) {
	if constexpr(!multi::has_dimensionality<std::decay_t<A>>::value) {
		return map(std::forward<F>(fun), identity_bind<A>{std::forward<A>(alpha)} ^ multi::extents_t<0>{}, std::forward<B>(omega));
	} else if constexpr(!multi::has_dimensionality<std::decay_t<B>>::value) {
		return map(std::forward<F>(fun), std::forward<A>(alpha), identity_bind<B>{std::forward<B>(omega)} ^ multi::extents_t<0>{});
	} else {
		using std::get;
		if constexpr(std::decay_t<A>::dimensionality < std::decay_t<B>::dimensionality) {
			return map(std::forward<F>(fun), std::forward<A>(alpha).repeated(get<std::decay_t<B>::dimensionality - std::decay_t<A>::dimensionality - 1>(omega.sizes())), std::forward<B>(omega));
		} else if constexpr(std::decay_t<B>::dimensionality < std::decay_t<A>::dimensionality) {
			return map(std::forward<F>(fun), std::forward<A>(alpha), std::forward<B>(omega).repeated(get<std::decay_t<A>::dimensionality - std::decay_t<B>::dimensionality - 1>(alpha.sizes())));
		} else {
			return apply(std::forward<F>(fun), std::forward<A>(alpha), std::forward<B>(omega));
		}
	}
}

// template<class A, class B>
// class apply_plus_t {
// 	A a_;  // NOLINT(cppcoreguidelines-avoid-const-or-ref-data-members)
// 	B b_;  // NOLINT(cppcoreguidelines-avoid-const-or-ref-data-members)

//  public:
// 	template<class AA, class BB>
// 	apply_plus_t(AA&& a, BB&& b) : a_{std::forward<AA>(a)}, b_{std::forward<BB>(b)} {}

// 	template<class... Is>
// 	constexpr auto operator()(Is... is) const {
// 		return multi::detail::invoke_square(a_, is...)  // like a_[is...] in C++23
// 			 +
// 			   multi::detail::invoke_square(b_, is...)  // like b_[is...] in C++23
// 			;
// 	}
// };

namespace detail {
struct plus {
	template<class T1, class T2>
	constexpr auto operator()(T1&& a, T2&& b) const;
};
}  // end namespace detail

template<class A, class B>
constexpr auto operator+(A&& alpha, B&& omega) /*noexcept*/ {
	return elementwise::map(elementwise::detail::plus{}, std::forward<A>(alpha), std::forward<B>(omega));
}

template<class T1, class T2>
constexpr auto detail::plus::operator()(T1&& a, T2&& b) const {
	using elementwise::operator+;  // cppcheck-suppress constStatement ;
	return std::forward<T1>(a) + std::forward<T2>(b);
}

namespace detail {
struct minus {
	template<class T1, class T2>
	constexpr auto operator()(T1&& a, T2&& b) const;
};
}  // end namespace detail

template<class A, class B>
constexpr auto operator-(A&& alpha, B&& omega) noexcept {
	return elementwise::map(elementwise::detail::minus{}, std::forward<A>(alpha), std::forward<B>(omega));
}

template<class T1, class T2>
constexpr auto detail::minus::operator()(T1&& a, T2&& b) const {
	using elementwise::operator-;  // cppcheck-suppress constStatement ;
	return std::forward<T1>(a) - std::forward<T2>(b);
}

template<class A>
constexpr auto operator-(A&& alpha) { return elementwise::apply(std::negate<>{}, std::forward<A>(alpha)); }

template<class A, class B>
constexpr auto operator*(A&& alpha, B&& omega) { return elementwise::map(std::multiplies<>{}, std::forward<A>(alpha), std::forward<B>(omega)); }

template<class A, class B>
constexpr auto operator/(A&& alpha, B&& omega) {
	return elementwise::map(std::divides<>{}, std::forward<A>(alpha), std::forward<B>(omega));
}

template<class T = void>
struct default_zero_f {
	template<class TT = T>
	auto operator()(TT const& /*unused*/) const { return TT{}; }
};

template<class T, class ZF>
constexpr auto eye(multi::ssize_t size, T unit, ZF zero_f) {
	return restricted([unit, zero = zero_f(unit)](auto ii, auto jj) { return ii == jj ? unit : zero; }, multi::extents_t<2>({size, size}));
}

template<class T, class ZF = default_zero_f<T>>
constexpr auto eye(multi::ssize_t size, T unit) {
	return eye(size, unit, default_zero_f<T>{});
}

template<class T = int>
constexpr auto eye(multi::ssize_t size) {
	return eye(size, T{1});
}

template<class Array, class DefaultZero>
constexpr auto zeros(Array&& arr, DefaultZero df) {
	auto exts = arr.extents();  // mull-ignore: cxx_init_const
	return restricted([arr_ = std::forward<Array>(arr), df](auto... ijk) { return df(arr_(ijk...)); }, exts);
}

template<typename Element = int, class Array, class DefaultZero = default_zero_f<Element>>
constexpr auto zeros(Array&& arr) {
	return zeros(std::forward<Array>(arr), DefaultZero{});
}

template<typename Element, dimensionality_type D>
constexpr auto zeros(multi::extents_t<D> const& exts) {
	return zeros<Element, multi::extents_t<D> const&>(exts);
}

template<dimensionality_type D>
constexpr auto zeros(multi::extents_t<D> const& exts) {
	return zeros<int, D>(exts);
}

// constrained to multi expressions (operands with dimensionality) so this does not become a
// greedy ADL `operator&&` candidate for unrelated types whose template args pull in namespace multi
template<class A, class B, std::enable_if_t<has_dimensionality<std::decay_t<A>>::value || has_dimensionality<std::decay_t<B>>::value, int> = 0>  // NOLINT(modernize-use-constraints) TODO(correaa)
constexpr auto operator&&(A&& alpha, B&& omega) { return elementwise::map(std::logical_and<>{}, std::forward<A>(alpha), std::forward<B>(omega)); }

template<class A, class B, std::enable_if_t<has_dimensionality<std::decay_t<A>>::value || has_dimensionality<std::decay_t<B>>::value, int> = 0>  // NOLINT(modernize-use-constraints) TODO(correaa)
constexpr auto operator||(A&& alpha, B&& omega) { return elementwise::map(std::logical_or<>{}, std::forward<A>(alpha), std::forward<B>(omega)); }

template<class F, class A, std::enable_if_t<true, decltype(std::declval<F&&>()(std::declval<typename std::decay_t<A>::element>()))*> = nullptr>  // NOLINT(modernize-use-constraints) for C++23
constexpr auto operator|(A&& alpha, F fun) {
	return std::forward<A>(alpha).element_transformed(fun);
}

template<class F, class A, std::enable_if_t<true, decltype(std::declval<F>()(std::declval<typename std::decay_t<A>::reference>()))*> = nullptr>  // NOLINT(modernize-use-constraints) for C++23
constexpr auto operator|(A&& a, F fun) {
	return std::forward<A>(a).transformed(fun);
}

template<class A>
class exp_bind_t {
	A a_;  // NOLINT(cppcoreguidelines-avoid-const-or-ref-data-members) TODO(correaa) consider saving .home() cursor

 public:
	template<class AA, std::enable_if_t<!std::is_base_of_v<exp_bind_t, std::decay_t<AA>>, int> = 0>  // NOLINT(modernize-use-constraints) for C++20
	BOOST_MULTI_HD constexpr explicit exp_bind_t(AA&& a) noexcept : a_{std::forward<AA>(a)} {}       // NOLINT(bugprone-forwarding-reference-overload)

	template<class... Is>
	constexpr auto operator()(Is... is) const {
		using ::std::exp;
		return exp(multi::detail::invoke_square(a_, is...));  // a_[is...] in C++23
	}
};

template<class A> exp_bind_t(A) -> exp_bind_t<A>;

template<class A, std::enable_if_t<multi::has_extents<std::decay_t<A>>::value, int> = 0>  // NOLINT(modernize-use-constraints) for C++23
BOOST_MULTI_HD constexpr auto exp(A&& alpha) {
	// shouldn't get to this point for scalars
	auto xs = alpha.extents();  // mull-ignore: cxx_init_const
	return exp_bind_t<A>(std::forward<A>(alpha)) ^ xs;
}

template<class A>
class log_bind_t {
	A a_;  // NOLINT(cppcoreguidelines-avoid-const-or-ref-data-members) TODO(correaa) consider saving .home() cursor

 public:
	template<class AA, std::enable_if_t<!std::is_base_of_v<log_bind_t, std::decay_t<AA>>, int> = 0>  // NOLINT(modernize-use-constraints) for C++20
	BOOST_MULTI_HD constexpr explicit log_bind_t(AA&& a) noexcept : a_{std::forward<AA>(a)} {}       // NOLINT(bugprone-forwarding-reference-overload)

	template<class... Is>
	constexpr auto operator()(Is... is) const {
		using ::std::log;
		return log(multi::detail::invoke_square(a_, is...));  // a_[is...] in C++23
	}
};

template<class A> log_bind_t(A) -> log_bind_t<A>;

template<class A, std::enable_if_t<multi::has_extents<std::decay_t<A>>::value, int> = 0>  // NOLINT(modernize-use-constraints) for C++23
BOOST_MULTI_HD constexpr auto log(A&& alpha) {
	auto xs = alpha.extensions();  // shouldn't get to this point for scalars
	return log_bind_t<A>(std::forward<A>(alpha)) ^ xs;
}

template<class A>
struct abs_bind_t {
	A a_;  // NOLINT(cppcoreguidelines-avoid-const-or-ref-data-members)

	template<class... Is>
	constexpr auto operator()(Is... is) const {
		using ::std::abs;
		return abs(multi::detail::invoke_square(a_, is...));  // a_[is...] in C++23
	}
};

template<class A>
constexpr auto abs(A const& a) { return abs_bind_t<decltype(a.home())>{a.home()} ^ a.extents(); }

template<class T> constexpr auto abs(std::initializer_list<T> il) { return abs(multi::array<T, 1>{il}); }
template<class T> constexpr auto abs(std::initializer_list<std::initializer_list<T>> il) { return abs(multi::array<T, 2>{il}); }

// #endif
}  // end namespace elementwise

// namespace broadcast = elementwise;

}  // end namespace boost::multi

#undef BOOST_MULTI_HD

#endif  // BOOST_MULTI_ELEMENTWISE_HPP
