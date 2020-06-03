#ifdef COMPILATION_INSTRUCTIONS//-*-indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4;-*-
$CXX -std=c++14 -D_TEST_MULTI_COMPLEX -x c++ $0 -o $0x&&$0x&&rm $0x;exit
#endif
// © Alfredo Correa 2020
#ifndef MULTI_COMPLEX_HPP
#define MULTI_COMPLEX_HPP

#include "array_ref.hpp"

#include<complex>

namespace boost{
namespace multi{

struct real_t;
struct imag_t;

template<class ValueType = double>
struct complex{
	using value_type = ValueType;
	value_type real;
	value_type imag;
	complex() = default;
	complex(value_type real) : real{real}, imag{value_type{0}}{}
	complex(value_type real, value_type imag) : real{real}, imag{imag}{}
	complex(std::complex<ValueType> const& other) : real{other.real()}, imag{other.imag()}{}
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
	operator T const&() const&{return reinterpret_cast<T const&>(*this);}
	template<
		class T,
		std::enable_if_t<
			sizeof(T)==2*sizeof(value_type) and 
			std::is_assignable<typename T::value_type&, decltype(std::declval<T>().real())>{} and
			std::is_assignable<typename T::value_type&, decltype(std::declval<T>().imag())>{}, int
		> =0
	>
	operator T&()&{return reinterpret_cast<T const&>(*this);}
	std::complex<value_type> const& std() const&{return reinterpret_cast<std::complex<value_type> const&>(*this);}
	std::complex<value_type>& std()&{return reinterpret_cast<std::complex<value_type>&>(*this);}
	friend auto abs(complex const& self){return abs(self.std());}
	friend complex operator-(complex const& self, complex const& other){return self.std() - other.std();}
};

struct real_t{
	template<class Array, typename E = typename std::decay_t<Array>::element, typename ValueType = typename E::value_type> 
	auto operator()(Array&& a) const
	->decltype(multi::member_array_cast<ValueType>(multi::reinterpret_array_cast<complex<ValueType>>(std::forward<Array>(a)), &complex<ValueType>::real)){
		return multi::member_array_cast<ValueType>(multi::reinterpret_array_cast<complex<ValueType>>(std::forward<Array>(a)), &complex<ValueType>::real);}
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
	->decltype(multi::member_array_cast<ValueType>(multi::reinterpret_array_cast<complex<ValueType>>(std::forward<Array>(a)), &complex<ValueType>::imag)){
		return multi::member_array_cast<ValueType>(multi::reinterpret_array_cast<complex<ValueType>>(std::forward<Array>(a)), &complex<ValueType>::imag);}
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

static real_t const real __attribute__((unused)) ;
static imag_t const imag __attribute__((unused)) ;

}}

namespace std{
	template<class T>
	struct is_trivially_default_constructible<std::complex<T>> : 
		is_trivially_default_constructible<T>{};
}

#if _TEST_MULTI_COMPLEX

#include<cassert>
#include<complex>
#include "array.hpp"

namespace multi = boost::multi;

template<class T> void what(T&&)=delete;

int main(){
	static_assert( std::is_trivially_default_constructible<std::complex<double>>{}, "!");
	static_assert( std::is_trivially_copy_constructible<std::complex<double>>{}, "!");
	static_assert( std::is_trivially_assignable<std::complex<double>&, std::complex<double> const>{}, "!");

	using complex = std::complex<double>;

	multi::array<complex, 2> A = {
		{ {1.,2.}, {3.,4.} },
		{ {22.,33.}, {5.,9.} }
	};

	{
		auto&& Areal = multi::member_array_cast<double>(A, &multi::complex<double>::real);
		auto&& Aimag = multi::member_array_cast<double>(A, &multi::complex<double>::imag);

		assert( Areal[1][0] == 22. );
		assert( Aimag[1][0] == 33. );
	}
	{
		auto&& Areal = multi::member_array_cast<double>(A, &multi::complex<double>::real);
		auto&& Aimag = multi::member_array_cast<double>(A, &multi::complex<double>::imag);

		assert( Areal[1][0] == 22. );
		assert( Aimag[1][0] == 33. );
	}
	{
		auto&& Areal = multi::real(A); // multi::real(A);
		auto&& Aimag = multi::imag(A); // multi::real(A);

		assert( &Areal[1][0] == &multi::real(A[1][0]) );
		assert( &Aimag[1][0] == &multi::imag(A[1][0]) );
	}
}

#endif
#endif


