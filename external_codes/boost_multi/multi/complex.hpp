#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;-*-
$CXX $0 -o $0x&&$0x&&rm $0x;exit
#endif
// © Alfredo Correa 2020

#ifndef MULTI_COMPLEX_HPP
#define MULTI_COMPLEX_HPP

#include "array_ref.hpp"

#include<complex>

namespace boost{
namespace multi{

constexpr class adl_conj_fn__{
	template<class... As>          auto _(priority<1>,        As&&... as) const JUSTRETURN(              std::conj(std::forward<As>(as)...))
	template<class... As>          auto _(priority<2>,        As&&... as) const DECLRETURN(                   conj(std::forward<As>(as)...))
	template<class T, class... As> auto _(priority<3>, T&& t, As&&... as) const DECLRETURN(std::forward<T>(t).conj(std::forward<As>(as)...))
public:
	template<class... As> auto operator()(As&&... as) const DECLRETURN(_(priority<3>{}, std::forward<As>(as)...))
} adl_conj;

constexpr class adl_real_fn__{
	template<class... As>          auto _(priority<1>,        As&&... as) const DECLRETURN(              std::real(std::forward<As>(as)...))
	template<class... As>          auto _(priority<2>,        As&&... as) const DECLRETURN(                   real(std::forward<As>(as)...))
	template<class T, class... As> auto _(priority<3>, T&& t, As&&... as) const DECLRETURN(std::forward<T>(t).real(std::forward<As>(as)...))
public:
	template<class... As> auto operator()(As&&... as) const DECLRETURN(_(priority<3>{}, std::forward<As>(as)...))
} adl_real;

constexpr class adl_imag_fn__{
	template<class... As>          auto _(priority<1>,        As&&... as) const DECLRETURN(              std::imag(std::forward<As>(as)...))
	template<class... As>          auto _(priority<2>,        As&&... as) const DECLRETURN(                   imag(std::forward<As>(as)...))
	template<class T, class... As> auto _(priority<3>, T&& t, As&&... as) const DECLRETURN(std::forward<T>(t).imag(std::forward<As>(as)...))
public:
	template<class... As> auto operator()(As&&... as) const DECLRETURN(_(priority<3>{}, std::forward<As>(as)...))
} adl_imag;

struct real_t;
struct imag_t;

template<class ValueType = double>
struct complex{
	using value_type = ValueType;
	value_type re;
	value_type im;
	complex() = default;
	// cppcheck-suppress noExplicitConstructor ; a real is a special complex without loss
	constexpr complex(value_type real) : re{real}, im{value_type{0}}{}
	constexpr complex(value_type real, value_type imag) : re{real}, im{imag}{}
	// cppcheck-suppress noExplicitConstructor ; can be copied from standard type
	constexpr complex(std::complex<ValueType> const& other) : re{other.real()}, im{other.imag()}{}
/*	friend value_type const& real(complex const& c){return c.real;}
	friend value_type      & real(complex      & c){return c.real;}
	friend value_type const& imag(complex const& c){return c.imag;}
	friend value_type      & imag(complex      & c){return c.imag;}*/
	template<
		class T,
		std::enable_if_t<
			sizeof(T)==2*sizeof(value_type) and 
			std::is_assignable<typename T::value_type&, decltype(std::declval<T>().real())>{} and
			std::is_assignable<typename T::value_type&, decltype(std::declval<T>().imag())>{}, int
		> =0
	>
	constexpr operator T const&() const&{return reinterpret_cast<T const&>(*this);}
	template<
		class T,
		std::enable_if_t<
			sizeof(T)==2*sizeof(value_type) and 
			std::is_assignable<typename T::value_type&, decltype(std::declval<T>().real())>{} and
			std::is_assignable<typename T::value_type&, decltype(std::declval<T>().imag())>{}, int
		> =0
	>
	constexpr operator T&()&{return reinterpret_cast<T const&>(*this);}
	constexpr std::complex<value_type> const& std() const&{return reinterpret_cast<std::complex<value_type> const&>(*this);}
	constexpr std::complex<value_type>      & std()      &{return reinterpret_cast<std::complex<value_type>      &>(*this);}
	friend constexpr auto abs(complex const& self){return abs(self.std());}
	friend constexpr complex operator-(complex const& self, complex const& other){return self.std() - other.std();}
	constexpr value_type& real()&{return re;}
	constexpr value_type const& real() const&{return re;}
//	friend constexpr value_type& real(complex& self){return self.re;}
//	friend constexpr value_type const& real(complex const& self){return self.re;}

	constexpr value_type& imag()&{return im;}
	constexpr value_type const& imag() const&{return im;}
//	friend constexpr value_type& imag(complex& self){return self.im;}
//	friend constexpr value_type const& imag(complex const& self){return self.im;}
	
	template<class Real> constexpr auto operator+=(Real const& other)&->decltype(re += other, *this){return re += other, *this;}
	template<class Real> constexpr auto operator-=(Real const& other)&->decltype(re -= other, *this){return re -= other, *this;}
	template<class Real> constexpr auto operator*=(Real const& other)&->decltype(re *= other, im *= other, *this){return re *= other, im *= other, *this;}
	template<class Real> constexpr auto operator/=(Real const& other)&->decltype(re /= other, im /= other, *this){return re /= other, im /= other, *this;}

	template<class Complex> constexpr auto operator+=(Complex const& other)&->decltype(re += other.re, im += other.im, *this){return re += other.re, im += other.im, *this;}
	template<class Complex> constexpr auto operator-=(Complex const& other)&->decltype(re -= other.re, im -= other.im, *this){return re -= other.re, im -= other.im, *this;}
};

struct real_t{
	template<class Array, typename E = typename std::decay_t<Array>::element, typename ValueType = typename E::value_type> 
	auto operator()(Array&& a) const
	->decltype(std::forward<Array>(a).template reinterpret_array_cast<complex<ValueType>>().template member_cast<ValueType>(&complex<ValueType>::real)){
		return std::forward<Array>(a).template reinterpret_array_cast<complex<ValueType>>().template member_cast<ValueType>(&complex<ValueType>::real);}
	template<class T, typename ValueType = typename std::decay_t<T>::value_type,
		std::enable_if_t<
			sizeof(T)==2*sizeof(ValueType) and 
			std::is_assignable<ValueType&, decltype(real(std::declval<T>()))>{} and
			std::is_assignable<ValueType&, decltype(imag(std::declval<T>()))>{}, int
		> =0
	>
	ValueType& operator()(T& t) const{return reinterpret_cast<multi::complex<ValueType>&>(t).real;}
	template<class T, typename ValueType = typename std::decay_t<T>::value_type,
		std::enable_if_t<
			sizeof(T)==2*sizeof(ValueType) and 
			std::is_assignable<ValueType&, decltype(real(std::declval<T>()))>{} and
			std::is_assignable<ValueType&, decltype(imag(std::declval<T>()))>{}, int
		> =0
	>
	ValueType const& operator()(T const& t) const{return reinterpret_cast<multi::complex<ValueType> const&>(t).real;}
};

struct imag_t{
	template<class Array, typename E = typename std::decay_t<Array>::element, typename ValueType = typename E::value_type> 
	auto operator()(Array&& a) const
	->decltype(std::forward<Array>(a).template reinterpret_array_cast<complex<ValueType>>().template member_cast<ValueType>(&complex<ValueType>::imag)){
		return std::forward<Array>(a).template reinterpret_array_cast<complex<ValueType>>().template member_cast<ValueType>(&complex<ValueType>::imag);}
	template<class T, typename ValueType = typename std::decay_t<T>::value_type, 
		std::enable_if_t<
			sizeof(T)==2*sizeof(ValueType) and 
			std::is_assignable<ValueType&, decltype(real(std::declval<T>()))>{} and
			std::is_assignable<ValueType&, decltype(imag(std::declval<T>()))>{}, int
		> =0
	>
	ValueType& operator()(T& t) const{return reinterpret_cast<multi::complex<ValueType>&>(t).imag;}
	template<class T, typename ValueType = typename std::decay_t<T>::value_type,
		std::enable_if_t<
			sizeof(T)==2*sizeof(ValueType) and 
			std::is_assignable<ValueType&, decltype(real(std::declval<T>()))>{} and
			std::is_assignable<ValueType&, decltype(imag(std::declval<T>()))>{}, int
		> =0
	>
	ValueType const& operator()(T const& t) const{return reinterpret_cast<multi::complex<ValueType> const&>(t).imag;}
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

}}

namespace std{
	template<class T>
	struct is_trivially_default_constructible<std::complex<T>> : 
		is_trivially_default_constructible<T>{};
}

#if not __INCLUDE_LEVEL__ // _TEST_MULTI_COMPLEX

#include<cassert>
#include<complex>
#include "array.hpp"

namespace multi = boost::multi;

template<class T> void what(T&&)=delete;

int main(){
	static_assert( std::is_trivially_default_constructible<std::complex<double>>{}, "!");
	static_assert( std::is_trivially_copy_constructible<std::complex<double>>{}, "!");
	static_assert( std::is_trivially_assignable<std::complex<double>&, std::complex<double> const>{}, "!");

	using complex = multi::complex<double>;

	multi::array<complex, 2> A = {
		{ {1.,2.}, {3.,4.} },
		{ {22.,33.}, {5.,9.} }
	};

	{
		auto&& Areal = A.member_cast<double>(&multi::complex<double>::re);
		auto&& Aimag = A.member_cast<double>(&multi::complex<double>::im);

		assert( Areal[1][0] == 22. );
		assert( Aimag[1][0] == 33. );
	}
	{
		auto&& Areal = A.member_cast<double>(&multi::complex<double>::re);
		auto&& Aimag = A.member_cast<double>(&multi::complex<double>::im);

		assert( Areal[1][0] == 22. );
		assert( Aimag[1][0] == 33. );
	}
	{
//		auto&& Areal = multi::real(A); // multi::real(A);
//		auto&& Aimag = multi::imag(A); // multi::real(A);

//		assert( &Areal[1][0] == &multi::real(A[1][0]) );
//		assert( &Aimag[1][0] == &multi::imag(A[1][0]) );
	}
}

#endif
#endif


