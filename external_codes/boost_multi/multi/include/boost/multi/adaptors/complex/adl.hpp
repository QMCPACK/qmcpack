// Copyright 2020-2024 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// TODO(correaa) move this header to blas/numeric

#ifndef BOOST_MULTI_ADAPTORS_COMPLEX_ADL_HPP
#define BOOST_MULTI_ADAPTORS_COMPLEX_ADL_HPP
#pragma once

// #include <boost/multi/array_ref.hpp>

// #include "detail/fix_complex_traits.hpp"

#include <complex>
#include <utility>  // for forward

#define BOOST_MULTI_DECLRETURN(ExpR) -> decltype(ExpR) {return ExpR;}  // NOLINT(cppcoreguidelines-macro-usage) saves a lot of typing
#define BOOST_MULTI_JUSTRETURN(ExpR)                   {return ExpR;}  // NOLINT(cppcoreguidelines-macro-usage) saves a lot of typing

namespace boost {  // NOLINT(modernize-concat-nested-namespaces) keep c++14 compat
namespace multi {

class adl_conj_t {
	template<class... As> constexpr auto _(priority<1> /**/, As&&... args) const BOOST_MULTI_JUSTRETURN(std::conj(std::forward<As>(args)...)) template<class... As> constexpr auto _(priority<2> /**/, As&&... args) const BOOST_MULTI_DECLRETURN(conj(std::forward<As>(args)...)) template<class T, class... As> constexpr auto _(priority<3> /**/, T&& arg, As&&... args) const BOOST_MULTI_DECLRETURN(std::forward<T>(arg).conj(std::forward<As>(args)...))

		public : template<class... As>
		         constexpr auto
		         operator()(As&&... args) const BOOST_MULTI_DECLRETURN(_(priority<3>{}, std::forward<As>(args)...))
};

inline constexpr adl_conj_t adl_conj;

class adl_real_t {
	template<class... As> constexpr auto _(priority<1> /**/, As&&... args) const BOOST_MULTI_DECLRETURN(std::real(std::forward<As>(args)...)) template<class... As> constexpr auto _(priority<2> /**/, As&&... args) const BOOST_MULTI_DECLRETURN(real(std::forward<As>(args)...)) template<class T, class... As> constexpr auto _(priority<3> /**/, T&& arg, As&&... args) const BOOST_MULTI_DECLRETURN(std::forward<T>(arg).real(std::forward<As>(args)...))

		public : template<class... As>
		         constexpr auto
		         operator()(As&&... args) const BOOST_MULTI_DECLRETURN(_(priority<3>{}, std::forward<As>(args)...))
};

inline constexpr adl_real_t adl_real;

class adl_imag_t {
	template<class... As> constexpr auto _(priority<1> /**/, As&&... args) const BOOST_MULTI_DECLRETURN(std::imag(std::forward<As>(args)...)) template<class... As> constexpr auto _(priority<2> /**/, As&&... args) const BOOST_MULTI_DECLRETURN(imag(std::forward<As>(args)...)) template<class T, class... As> constexpr auto _(priority<3> /**/, T&& arg, As&&... args) const BOOST_MULTI_DECLRETURN(std::forward<T>(arg).imag(std::forward<As>(args)...))

		public : template<class... As>
		         constexpr auto
		         operator()(As&&... args) const BOOST_MULTI_DECLRETURN(_(priority<3>{}, std::forward<As>(args)...))
};

inline constexpr adl_imag_t adl_imag;

struct real_t;
struct imag_t;

template<class ValueType = double>
struct _complex {  // NOLINT(readability-identifier-naming) deprecating this
	using value_type = ValueType;

 private:
	value_type re_;
	value_type im_;

 public:
	_complex() = default;

	constexpr explicit _complex(value_type real) : re_{real}, im_{value_type{0}} {}
	constexpr _complex(value_type real, value_type imag)  // NOLINT(bugprone-easily-swappable-parameters)
	: re_{real}, im_{imag} {}

	constexpr explicit _complex(std::complex<ValueType> const& other) : re_{other.real()}, im_{other.imag()} {}

	template<
		class T,
		std::enable_if_t<  // NOLINT(modernize-use-constraints) TODO(correaa) for C++20
			sizeof(T) == 2 * sizeof(value_type) &&
				std::is_assignable<typename T::value_type&, decltype(std::declval<T>().real())>{} && std::is_assignable<typename T::value_type&, decltype(std::declval<T>().imag())>{},
			int> =0>
	constexpr explicit operator T const&() const& {
		return reinterpret_cast<T const&>(*this);  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast)
	}
	template<
		class T,
		std::enable_if_t<  // NOLINT(modernize-use-constraints) TODO(correaa) for C++20
			sizeof(T) == 2 * sizeof(value_type) &&
				std::is_assignable_v<typename T::value_type&, decltype(std::declval<T>().real())> &&
				std::is_assignable_v<typename T::value_type&, decltype(std::declval<T>().imag())>,
			int> = 0>
	constexpr explicit operator T&() & { return reinterpret_cast<T const&>(*this); }  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast)

	constexpr auto std() const& -> std::complex<value_type> const& {
		return reinterpret_cast<std::complex<value_type> const&>(*this);  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast)
	}
	constexpr auto std() & -> std::complex<value_type>& {
		return reinterpret_cast<std::complex<value_type>&>(*this);  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast)
	}

	friend constexpr auto abs(_complex const& self) { return abs(self.std()); }
	friend constexpr auto operator-(_complex const& self, _complex const& other)
		-> _complex { return self.std() - other.std(); }

	constexpr auto real() & -> value_type& { return re_; }
	constexpr auto real() const& -> value_type const& { return re_; }

	constexpr auto imag() & -> value_type& { return im_; }
	constexpr auto imag() const& -> value_type const& { return im_; }

	template<class Real> constexpr auto operator+=(Real const& other) & -> decltype(re_ += other, *this) { return re_ += other, *this; }
	template<class Real> constexpr auto operator-=(Real const& other) & -> decltype(re_ -= other, *this) { return re_ -= other, *this; }
	template<class Real> constexpr auto operator*=(Real const& other) & -> decltype(re_ *= other, im_ *= other, *this) { return re_ *= other, im_ *= other, *this; }
	template<class Real> constexpr auto operator/=(Real const& other) & -> decltype(re_ /= other, im_ /= other, *this) { return re_ /= other, im_ /= other, *this; }

	template<class Complex> constexpr auto operator+=(Complex const& other) & -> decltype(re_ += other.re, im_ += other.im, *this) { return re_ += other.re, im_ += other.im, *this; }
	template<class Complex> constexpr auto operator-=(Complex const& other) & -> decltype(re_ -= other.re, im_ -= other.im, *this) { return re_ -= other.re, im_ -= other.im, *this; }
};

struct real_t {
	template<class Array, typename E = typename std::decay_t<Array>::element, typename ValueType = typename E::value_type>
	constexpr auto operator()(Array&& array) const
		-> decltype(std::forward<Array>(array).template reinterpret_array_cast<_complex<ValueType>>().template member_cast<ValueType>(&_complex<ValueType>::real)) {
		return std::forward<Array>(array).template reinterpret_array_cast<_complex<ValueType>>().template member_cast<ValueType>(&_complex<ValueType>::real);
	}
	template<class T, typename ValueType = typename std::decay_t<T>::value_type,
	         std::enable_if_t<  // NOLINT(modernize-use-constraints) TODO(correaa) for C++20
		         sizeof(T) == 2 * sizeof(ValueType) &&
			         std::is_assignable_v<ValueType&, decltype(real(std::declval<T>()))> &&
			         std::is_assignable_v<ValueType&, decltype(imag(std::declval<T>()))>,
		         int> = 0>
	constexpr auto operator()(T& value) const -> ValueType& { return reinterpret_cast<multi::_complex<ValueType>&>(value).real; }  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast) : TODO(correaa) : t[0]
	template<class T, typename ValueType = typename std::decay_t<T>::value_type,
	         std::enable_if_t<  // NOLINT(modernize-use-constraints) TODO(correaa) for C++20
		         sizeof(T) == 2 * sizeof(ValueType) &&
			         std::is_assignable_v<ValueType&, decltype(real(std::declval<T>()))> &&
			         std::is_assignable_v<ValueType&, decltype(imag(std::declval<T>()))>,
		         int> = 0>
	auto operator()(T const& value) const -> ValueType const& {
		return reinterpret_cast<multi::_complex<ValueType> const&>(value).real;  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast) : TODO(correaa) : t[0]
	}
};

struct imag_t {
	template<class Array, typename E = typename std::decay_t<Array>::element, typename ValueType = typename E::value_type>
	constexpr auto operator()(Array&& array) const
		-> decltype(std::forward<Array>(array).template reinterpret_array_cast<_complex<ValueType>>().template member_cast<ValueType>(&_complex<ValueType>::imag)) {
		return std::forward<Array>(array).template reinterpret_array_cast<_complex<ValueType>>().template member_cast<ValueType>(&_complex<ValueType>::imag);
	}
	template<class T, typename ValueType = typename std::decay_t<T>::value_type,
	         std::enable_if_t<  // NOLINT(modernize-use-constraints) TODO(correaa) for C++20
		         sizeof(T) == 2 * sizeof(ValueType) &&
			         std::is_assignable<ValueType&, decltype(real(std::declval<T>()))>{} &&
			         std::is_assignable<ValueType&, decltype(imag(std::declval<T>()))>{},
		         int> = 0>
	constexpr auto operator()(T& value) const -> ValueType& {
		return reinterpret_cast<multi::_complex<ValueType>&>(value).imag;  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast) : TODO(correaa) : t[1]
	}
	template<class T, typename ValueType = typename std::decay_t<T>::value_type,
	         std::enable_if_t<  // NOLINT(modernize-use-constraints) TODO(correaa) for C++20
		         sizeof(T) == 2 * sizeof(ValueType) &&
			         std::is_assignable<ValueType&, decltype(real(std::declval<T>()))>{} &&
			         std::is_assignable<ValueType&, decltype(imag(std::declval<T>()))>{},
		         int> = 0>
	constexpr auto operator()(T const& value) const -> ValueType const& {
		return reinterpret_cast<multi::_complex<ValueType> const&>(value).imag;  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast) : TODO(correaa) : t[1]
	}
};

//[[maybe_unused]] static constexpr real_t real;
//[[maybe_unused]] static constexpr imag_t imag;

}  // end namespace multi
}  // end namespace boost

#undef BOOST_MULTI_DECLRETURN
#undef BOOST_MULTI_JUSTRETURN

#endif  // BOOST_MULTI_ADAPTORS_COMPLEX_ADL_HPP
