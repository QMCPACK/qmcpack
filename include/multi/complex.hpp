// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2020-2022 Alfredo Correa

#ifndef MULTI_COMPLEX_HPP
#define MULTI_COMPLEX_HPP

#include "array_ref.hpp"

#include "detail/fix_complex_traits.hpp"

#include<complex>
#include<utility>  // for forward

namespace boost {  // NOLINT(modernize-concat-nested-namespaces) keep c++14 compat
namespace multi {

constexpr class adl_conj_t {
	template<class... As>          constexpr auto _(priority<1>/**/,          As&&... args) const JUSTRETURN(              std::  conj(std::forward<As>(args)...))
	template<class... As>          constexpr auto _(priority<2>/**/,          As&&... args) const DECLRETURN(                     conj(std::forward<As>(args)...))
	template<class T, class... As> constexpr auto _(priority<3>/**/, T&& arg, As&&... args) const DECLRETURN(std::forward<T>(arg).conj(std::forward<As>(args)...))

 public:
	template<class... As> constexpr auto operator()(As&&... args) const DECLRETURN(_(priority<3>{}, std::forward<As>(args)...))
} adl_conj;

constexpr class adl_real_t {
	template<class... As>          constexpr auto _(priority<1>/**/,          As&&... args) const DECLRETURN(                std::real(std::forward<As>(args)...))
	template<class... As>          constexpr auto _(priority<2>/**/,          As&&... args) const DECLRETURN(                     real(std::forward<As>(args)...))
	template<class T, class... As> constexpr auto _(priority<3>/**/, T&& arg, As&&... args) const DECLRETURN(std::forward<T>(arg).real(std::forward<As>(args)...))

 public:
	template<class... As> constexpr auto operator()(As&&... args) const DECLRETURN(_(priority<3>{}, std::forward<As>(args)...))
} adl_real;

constexpr class adl_imag_t {
	template<class... As>          constexpr auto _(priority<1>/**/,          As&&... args) const DECLRETURN(                std::imag(std::forward<As>(args)...))
	template<class... As>          constexpr auto _(priority<2>/**/,          As&&... args) const DECLRETURN(                     imag(std::forward<As>(args)...))
	template<class T, class... As> constexpr auto _(priority<3>/**/, T&& arg, As&&... args) const DECLRETURN(std::forward<T>(arg).imag(std::forward<As>(args)...))

 public:
	template<class... As> constexpr auto operator()(As&&... args) const DECLRETURN(_(priority<3>{}, std::forward<As>(args)...))
} adl_imag;

struct real_t;
struct imag_t;

template<class ValueType = double>
struct complex {
	using value_type = ValueType;

 private:
	value_type re;
	value_type im;

 public:
	complex() = default;

	constexpr explicit complex(value_type real) : re{real}, im{value_type{0}} {}
	constexpr complex(value_type real, value_type imag) // NOLINT(bugprone-easily-swappable-parameters)
	: re{real}, im{imag} {}

	constexpr explicit complex(std::complex<ValueType> const& other) : re{other.real()}, im{other.imag()} {}

	template<
		class T,
		std::enable_if_t<
			sizeof(T)==2*sizeof(value_type) and
			std::is_assignable<typename T::value_type&, decltype(std::declval<T>().real())>{} and
			std::is_assignable<typename T::value_type&, decltype(std::declval<T>().imag())>{}, int
		> =0
	>
	constexpr explicit operator T const&() const& {
		return reinterpret_cast<T const&>(*this);  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast)
	}
	template<
		class T,
		std::enable_if_t<
			sizeof(T)==2*sizeof(value_type) and
			std::is_assignable<typename T::value_type&, decltype(std::declval<T>().real())>{} and
			std::is_assignable<typename T::value_type&, decltype(std::declval<T>().imag())>{}, int
		> = 0
	>
	constexpr explicit operator T&()& {return reinterpret_cast<T const&>(*this);}  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast)

	constexpr auto std() const& -> std::complex<value_type> const& {
		return reinterpret_cast<std::complex<value_type> const&>(*this);  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast)
	}
	constexpr auto std()      & -> std::complex<value_type>      & {
		return reinterpret_cast<std::complex<value_type>      &>(*this);  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast)
	}

	friend constexpr auto abs(complex const& self) {return abs(self.std());}
	friend constexpr auto operator-(complex const& self, complex const& other)
	-> complex{return self.std() - other.std();}

	constexpr auto real()      & -> value_type      & {return re;}
	constexpr auto real() const& -> value_type const& {return re;}

	constexpr auto imag()      & -> value_type      & {return im;}
	constexpr auto imag() const& -> value_type const& {return im;}

	template<class Real> constexpr auto operator+=(Real const& other)&->decltype(re += other, *this) {return re += other, *this;}
	template<class Real> constexpr auto operator-=(Real const& other)&->decltype(re -= other, *this) {return re -= other, *this;}
	template<class Real> constexpr auto operator*=(Real const& other)&->decltype(re *= other, im *= other, *this) {return re *= other, im *= other, *this;}
	template<class Real> constexpr auto operator/=(Real const& other)&->decltype(re /= other, im /= other, *this) {return re /= other, im /= other, *this;}

	template<class Complex> constexpr auto operator+=(Complex const& other)&->decltype(re += other.re, im += other.im, *this) {return re += other.re, im += other.im, *this;}
	template<class Complex> constexpr auto operator-=(Complex const& other)&->decltype(re -= other.re, im -= other.im, *this) {return re -= other.re, im -= other.im, *this;}
};

struct real_t {
	template<class Array, typename E = typename std::decay_t<Array>::element, typename ValueType = typename E::value_type>
	constexpr auto operator()(Array&& array) const
	->decltype(std::forward<Array>(array).template reinterpret_array_cast<complex<ValueType>>().template member_cast<ValueType>(&complex<ValueType>::real)) {
		return std::forward<Array>(array).template reinterpret_array_cast<complex<ValueType>>().template member_cast<ValueType>(&complex<ValueType>::real); }
	template<class T, typename ValueType = typename std::decay_t<T>::value_type,
		std::enable_if_t<
			sizeof(T)==2*sizeof(ValueType) and
			std::is_assignable<ValueType&, decltype(real(std::declval<T>()))>{} and
			std::is_assignable<ValueType&, decltype(imag(std::declval<T>()))>{}, int
		> =0
	>
	constexpr auto operator()(T& value) const -> ValueType& {return reinterpret_cast<multi::complex<ValueType>&>(value).real;}  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast) : TODO(correaa) : t[0]
	template<class T, typename ValueType = typename std::decay_t<T>::value_type,
		std::enable_if_t<
			sizeof(T)==2*sizeof(ValueType) and
			std::is_assignable<ValueType&, decltype(real(std::declval<T>()))>{} and
			std::is_assignable<ValueType&, decltype(imag(std::declval<T>()))>{}, int
		> =0
	>
	auto operator()(T const& value) const -> ValueType const& {
		return reinterpret_cast<multi::complex<ValueType> const&>(value).real;  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast) : TODO(correaa) : t[0]
	}
};

struct imag_t {
	template<class Array, typename E = typename std::decay_t<Array>::element, typename ValueType = typename E::value_type>
	constexpr auto operator()(Array&& array) const
	->decltype(std::forward<Array>(array).template reinterpret_array_cast<complex<ValueType>>().template member_cast<ValueType>(&complex<ValueType>::imag)) {
		return std::forward<Array>(array).template reinterpret_array_cast<complex<ValueType>>().template member_cast<ValueType>(&complex<ValueType>::imag); }
	template<class T, typename ValueType = typename std::decay_t<T>::value_type,
		std::enable_if_t<
			sizeof(T) == 2*sizeof(ValueType) and
			std::is_assignable<ValueType&, decltype(real(std::declval<T>()))>{} and
			std::is_assignable<ValueType&, decltype(imag(std::declval<T>()))>{}, int
		> =0
	>
	constexpr auto operator()(T& value) const -> ValueType& {
		return reinterpret_cast<multi::complex<ValueType>&>(value).imag;  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast) : TODO(correaa) : t[1]
	}
	template<class T, typename ValueType = typename std::decay_t<T>::value_type,
		std::enable_if_t<
			sizeof(T)==2*sizeof(ValueType) and
			std::is_assignable<ValueType&, decltype(real(std::declval<T>()))>{} and
			std::is_assignable<ValueType&, decltype(imag(std::declval<T>()))>{}, int
		> =0
	>
	constexpr auto operator()(T const& value) const -> ValueType const&{
		return reinterpret_cast<multi::complex<ValueType> const&>(value).imag;  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast) : TODO(correaa) : t[1]
	}
};

[[maybe_unused]] static constexpr real_t real;
[[maybe_unused]] static constexpr imag_t imag;

}  // end namespace multi
}  // end namespace boost

static_assert( boost::multi::is_trivially_default_constructible<std::complex<double>>::value );
static_assert( boost::multi::is_trivially_default_constructible<std::complex<float >>::value );

static_assert( boost::multi::is_trivial<std::complex<double>>::value );
static_assert( boost::multi::is_trivial<std::complex<float >>::value );


#if defined(__INCLUDE_LEVEL__) and not __INCLUDE_LEVEL__

#include<cassert>
#include "array.hpp"

namespace multi = boost::multi;

template<class T> void what(T&&)=delete;

int main() {

	using complex = multi::complex<double>;

	multi::array<complex, 2> A = {
		{ {1. ,  2.}, {3., 4.} },
		{ {22., 33.}, {5., 9.} }
	};

	{
		auto&& Areal = A.member_cast<double>(&multi::complex<double>::re);
		auto&& Aimag = A.member_cast<double>(&multi::complex<double>::im);

		assert( Areal[1][0] == 22. );
		assert( Aimag[1][0] == 33. );
	} {
		auto&& Areal = A.member_cast<double>(&multi::complex<double>::re);
		auto&& Aimag = A.member_cast<double>(&multi::complex<double>::im);

		assert( Areal[1][0] == 22. );
		assert( Aimag[1][0] == 33. );
	}
}

#endif
#endif
