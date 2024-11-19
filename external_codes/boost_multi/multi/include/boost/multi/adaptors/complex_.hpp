// Copyright 2023-2024 Alfredo A. Correa

#ifndef BOOST_MULTI_ADAPTORS_COMPLEX_HPP
#define BOOST_MULTI_ADAPTORS_COMPLEX_HPP
#pragma once

#include<complex>  // to define its traits

namespace boost::multi {

template<class T> struct [[nodiscard]] complex;

template<class T> struct [[nodiscard]] imaginary;

template<class U>
constexpr auto operator+(U real, imaginary<U> imag) -> complex<U>;

template<class T>
struct [[nodiscard]] imaginary {
	T _value;  // NOLINT(misc-non-private-member-variables-in-classes) I want the class to be an aggregate

	using value_type = T;

	//  constexpr explicit imaginary(T value) : value_{value} {}
	template<class U>
	friend constexpr auto operator+(U real, imaginary<U> imag) -> complex<U>;
	//  constexpr static imaginary i{T{1}};  // NOLINT(clang-diagnostic-error) "constexpr variable cannot have non-literal type"?
	friend constexpr auto operator*(T real, imaginary imag) {
		return imaginary{real * imag._value};
	}
	[[nodiscard]] constexpr auto operator*(imaginary other) const { return -_value * other._value; }
	[[nodiscard]] constexpr auto operator/(imaginary other) const { return _value / other._value; }
	[[nodiscard]] constexpr auto operator+(imaginary other) const { return imaginary{_value + other._value}; }
	[[nodiscard]] constexpr auto operator-(imaginary other) const { return imaginary{_value + other._value}; }

	[[nodiscard]] constexpr auto operator==(imaginary const& other) const { return _value == other._value; };
	[[nodiscard]] constexpr auto operator!=(imaginary const& other) const { return _value != other._value; };
};

template<>
struct [[nodiscard]] imaginary<void> {
	template<class T>
	friend constexpr auto operator*(T real, imaginary /*self*/) { return imaginary<T>{real}; }
	template<class T>
	[[nodiscard]] constexpr auto operator*(imaginary<T> other) const { return -other._value; }
	template<class T>
	[[nodiscard]] constexpr auto operator/(imaginary<T> other) const { return T{1} / other._value; }
};

inline constexpr imaginary<void> I{};  // NOLINT(readability-identifier-length) imaginary unit

namespace literals {
// constexpr imaginary<double> operator""_i(unsigned long long d) {
//  return imaginary<double>{static_cast<double>(d)};
// }

constexpr auto operator""_i(long double value) { return imaginary<double>{static_cast<double>(value)}; }
//  constexpr auto operator""   i(long double value) {return imaginary<double>{static_cast<double>(value)};}
constexpr auto operator""_I(long double value) { return imaginary<double>{static_cast<double>(value)}; }

//  constexpr auto operator"" f_i(long double value) {return imaginary<float >{static_cast<float >(value)};}
constexpr auto operator""_f_i(long double value) { return imaginary<float>{static_cast<float>(value)}; }
constexpr auto operator""_if(long double value) { return imaginary<float>{static_cast<float>(value)}; }
constexpr auto operator""_F_I(long double value) { return imaginary<float>{static_cast<float>(value)}; }
constexpr auto operator""_IF(long double value) { return imaginary<float>{static_cast<float>(value)}; }

// template<char... Chars>
// constexpr auto operator""_FI() noexcept {}

}  // namespace literals

template<class T>
struct [[nodiscard]] complex {
	using real_type = T;

//  using value_type /*[[deprecated("reason")]]*/ = T;

	real_type _real;  // NOLINT(misc-non-private-member-variables-in-classes) complex should be an aggregate class
	real_type _imag;  // NOLINT(misc-non-private-member-variables-in-classes) complex should be an aggregate class

	template<class U>
	friend constexpr auto operator+(U real, imaginary<U> imag) -> complex<U>;

	friend constexpr auto operator*(real_type scale, complex self) { return complex{scale * self._real, scale * self._imag}; }
	friend constexpr auto operator/(complex self, real_type scale) { return complex{self._real / scale, self._imag / scale}; }

	friend constexpr auto operator+(T real, complex self) { return complex{real + self._real, self._imag}; }
	friend constexpr auto operator-(T real, complex self) { return complex{real - self._real, self._imag}; }

	[[nodiscard]] constexpr auto real() const -> real_type { return _real; }
	[[nodiscard]] constexpr auto imag() const -> real_type { return _imag; }

	friend constexpr auto conj(complex self) { return complex{self._real, -self._imag}; }

	constexpr auto operator==(complex const& other) const {return _real == other._real && _imag == other._imag;}
	constexpr auto operator!=(complex const& other) const {return _real != other._real || _imag != other._imag;}

//  auto operator=(complex const&) -> complex& = default;
	constexpr auto operator=(real_type re) -> complex& {(*this) = complex{re, real_type{0.0}}; return *this;}

	friend constexpr auto operator-(complex self) {return complex{-self._real, -self._imag};}
	friend constexpr auto operator+(complex self) {return complex{+self._real, +self._imag};}

	friend constexpr auto operator+(complex z1, complex z2) {return complex{z1._real + z2._real, z1._imag + z2._imag};}
	friend constexpr auto operator-(complex z1, complex z2) {return complex{z1._real - z2._real, z1._imag - z2._imag};}

	constexpr auto operator+=(complex other) -> complex& {_real += other._real; _imag += other._imag; return *this;}
	constexpr auto operator-=(complex other) -> complex& {_real -= other._real; _imag -= other._imag; return *this;}

	friend constexpr auto operator*(complex z1, complex z2) {
		return complex{z1._real * z2._real - z1._imag * z2._imag, z1._real * z2._imag + z1._imag * z2._real};
	}
	friend constexpr auto operator/(complex z1, complex z2) {
		auto const nrm = norm(z2);
		return complex{
			(z1._real * z2._real + z1._imag * z2._imag)/nrm,
			(z1._imag * z2._real - z1._real * z2._imag)/nrm
		};

		// typedef typename detail::promoted_numerical_type<T0, T1>::type T;

		// // Find `abs` by ADL.
		// using std::abs;

		// T s = abs(y.real()) + abs(y.imag());

		// T oos = T(1.0) / s;

		// T ars = x.real() * oos;
		// T ais = x.imag() * oos;
		// T brs = y.real() * oos;
		// T bis = y.imag() * oos;

		// s = (brs * brs) + (bis * bis);

		// oos = T(1.0) / s;

		// complex<T> quot( ((ars * brs) + (ais * bis)) * oos
		//              , ((ais * brs) - (ars * bis)) * oos);
		// return quot;
	}
	friend constexpr auto norm(complex self) {
		return self._real*self._real + self._imag*self._imag;  // TODO(correaa) revise this, use more exact formula
	}
	friend constexpr auto abs(complex self) {
	//  return hypot(z.real(), z.imag());
		using std::sqrt;
		return sqrt(self._real*self._real + self._real*self._real);  // bad! according to NR
		// using std::abs;
		// return self._real > self._imag?
		//       abs(self._real)*sqrt(real_type{1} + (self._imag/self._real)*(self._imag/self._real))
		//      :abs(self._imag)*sqrt(real_type{1} + (self._real/self._imag)*(self._real/self._imag))
		// ;
	}
};

template<class U>
constexpr auto operator+(U real, imaginary<U> imag) -> complex<U> { return {real, imag._value}; }

}  // end namespace boost::multi
#endif
