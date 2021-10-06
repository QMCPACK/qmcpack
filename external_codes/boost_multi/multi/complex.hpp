// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2020-2021 Alfredo Correa

#ifndef MULTI_COMPLEX_HPP
#define MULTI_COMPLEX_HPP

#include "array_ref.hpp"

#include<complex>
#include<utility>  // for forward

namespace boost {
namespace multi {

constexpr class adl_conj_t {
	template<class... As>          auto _(priority<1>/**/,        As&&... as) const JUSTRETURN(              std::conj(std::forward<As>(as)...))
	template<class... As>          auto _(priority<2>/**/,        As&&... as) const DECLRETURN(                   conj(std::forward<As>(as)...))
	template<class T, class... As> auto _(priority<3>/**/, T&& t, As&&... as) const DECLRETURN(std::forward<T>(t).conj(std::forward<As>(as)...))
public:
	template<class... As> auto operator()(As&&... as) const DECLRETURN(_(priority<3>{}, std::forward<As>(as)...))
} adl_conj;

constexpr class adl_real_t {
	template<class... As>          auto _(priority<1>/**/,        As&&... as) const DECLRETURN(              std::real(std::forward<As>(as)...))
	template<class... As>          auto _(priority<2>/**/,        As&&... as) const DECLRETURN(                   real(std::forward<As>(as)...))
	template<class T, class... As> auto _(priority<3>/**/, T&& t, As&&... as) const DECLRETURN(std::forward<T>(t).real(std::forward<As>(as)...))
public:
	template<class... As> auto operator()(As&&... as) const DECLRETURN(_(priority<3>{}, std::forward<As>(as)...))
} adl_real;

constexpr class adl_imag_t{
	template<class... As>          auto _(priority<1>/**/,        As&&... as) const DECLRETURN(              std::imag(std::forward<As>(as)...))
	template<class... As>          auto _(priority<2>/**/,        As&&... as) const DECLRETURN(                   imag(std::forward<As>(as)...))
	template<class T, class... As> auto _(priority<3>/**/, T&& t, As&&... as) const DECLRETURN(std::forward<T>(t).imag(std::forward<As>(as)...))
public:
	template<class... As> auto operator()(As&&... as) const DECLRETURN(_(priority<3>{}, std::forward<As>(as)...))
} adl_imag;

struct real_t;
struct imag_t;

template<class ValueType = double>
struct complex{
	using value_type = ValueType;

 private:
	value_type re;
	value_type im;

 public:
	complex() = default;

	constexpr explicit complex(value_type real) : re{real}, im{value_type{0}} {}
	constexpr complex(value_type real, value_type imag) : re{real}, im{imag} {}

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

struct real_t{
	template<class Array, typename E = typename std::decay_t<Array>::element, typename ValueType = typename E::value_type>
	auto operator()(Array&& a) const
	->decltype(std::forward<Array>(a).template reinterpret_array_cast<complex<ValueType>>().template member_cast<ValueType>(&complex<ValueType>::real)) {
		return std::forward<Array>(a).template reinterpret_array_cast<complex<ValueType>>().template member_cast<ValueType>(&complex<ValueType>::real); }
	template<class T, typename ValueType = typename std::decay_t<T>::value_type,
		std::enable_if_t<
			sizeof(T)==2*sizeof(ValueType) and
			std::is_assignable<ValueType&, decltype(real(std::declval<T>()))>{} and
			std::is_assignable<ValueType&, decltype(imag(std::declval<T>()))>{}, int
		> =0
	>
	auto operator()(T& t) const -> ValueType& {return reinterpret_cast<multi::complex<ValueType>&>(t).real;}  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast) : TODO(correaa) : t[0]
	template<class T, typename ValueType = typename std::decay_t<T>::value_type,
		std::enable_if_t<
			sizeof(T)==2*sizeof(ValueType) and
			std::is_assignable<ValueType&, decltype(real(std::declval<T>()))>{} and
			std::is_assignable<ValueType&, decltype(imag(std::declval<T>()))>{}, int
		> =0
	>
	auto operator()(T const& t) const -> ValueType const&{
		return reinterpret_cast<multi::complex<ValueType> const&>(t).real;  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast) : TODO(correaa) : t[0]
	}
};

struct imag_t {
	template<class Array, typename E = typename std::decay_t<Array>::element, typename ValueType = typename E::value_type>
	auto operator()(Array&& a) const
	->decltype(std::forward<Array>(a).template reinterpret_array_cast<complex<ValueType>>().template member_cast<ValueType>(&complex<ValueType>::imag)) {
		return std::forward<Array>(a).template reinterpret_array_cast<complex<ValueType>>().template member_cast<ValueType>(&complex<ValueType>::imag); }
	template<class T, typename ValueType = typename std::decay_t<T>::value_type,
		std::enable_if_t<
			sizeof(T)==2*sizeof(ValueType) and
			std::is_assignable<ValueType&, decltype(real(std::declval<T>()))>{} and
			std::is_assignable<ValueType&, decltype(imag(std::declval<T>()))>{}, int
		> = 0
	>
	auto operator()(T& t) const -> ValueType& {
		return reinterpret_cast<multi::complex<ValueType>&>(t).imag;  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast) : TODO(correaa) : t[1]
	}
	template<class T, typename ValueType = typename std::decay_t<T>::value_type,
		std::enable_if_t<
			sizeof(T)==2*sizeof(ValueType) and
			std::is_assignable<ValueType&, decltype(real(std::declval<T>()))>{} and
			std::is_assignable<ValueType&, decltype(imag(std::declval<T>()))>{}, int
		> =0
	>
	auto operator()(T const& t) const -> ValueType const&{
		return reinterpret_cast<multi::complex<ValueType> const&>(t).imag;  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast) : TODO(correaa) : t[1]
	}
};

/*
template<class T2, class Array, class P2 = typename std::pointer_traits<typename std::decay<Array>::type::element_ptr>::template rebind<T2>,
typename E = typename std::decay_t<Array>::element, 
typename R = decltype(real(E{})),
std::enable_if_t<sizeof(E)==2*sizeof(typename E::value_type), int> =0
>
decltype(auto) member_array_cast(Array&& a, real_t const*){
	struct Complex{double real; double imag;};
	return multi::member_array_cast<double>(multi::reinterpret_array_cast<Complex>(std::forward<Array>(a)), &Complex::real);
}
template<class T2, class Array, class P2 = typename std::pointer_traits<typename std::decay<Array>::type::element_ptr>::template rebind<T2>,
typename E = typename std::decay_t<Array>::element, 
typename R = decltype(real(E{})),
std::enable_if_t<sizeof(E)==2*sizeof(typename E::value_type), int> =0
>
decltype(auto) member_array_cast(Array&& a, imag_t const*){
	struct Complex{double real; double imag;};
	return multi::member_array_cast<double>(multi::reinterpret_array_cast<Complex>(std::forward<Array>(a)), &Complex::imag);
}
*/

static constexpr real_t real MAYBE_UNUSED;
static constexpr imag_t imag MAYBE_UNUSED;

}  // end namespace multi
}  // end namespace boost

namespace std {

template<class T>
struct is_trivially_default_constructible<std::complex<T>>
: is_trivially_default_constructible<T>{};

}  // end namespace std

#if defined(__INCLUDE_LEVEL__) and not __INCLUDE_LEVEL__

#include<cassert>c
#include "array.hpp"

namespace multi = boost::multi;

template<class T> void what(T&&)=delete;

int main() {
	static_assert( std::is_trivially_default_constructible<std::complex<double>>{}, "!");
	static_assert( std::is_trivially_copy_constructible<std::complex<double>>{}, "!");
	static_assert( std::is_trivially_assignable<std::complex<double>&, std::complex<double> const>{}, "!");

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


