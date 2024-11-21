// Copyright 2019-2024 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_MULTI_ADAPTORS_BLAS_OPERATIONS_HPP
#define BOOST_MULTI_ADAPTORS_BLAS_OPERATIONS_HPP

#include <boost/multi/adaptors/blas/numeric.hpp>

namespace boost::multi::blas {

// template<class M> auto transposed(M&& arr) -> decltype(auto) {return rotated(std::forward<M>(arr));}

template<class A>
[[deprecated("use blas::T to avoid conflict with multi::transposed")]]
auto transposed(A&& arr) -> decltype(auto) {return std::forward<A>(arr).rotated();}

template<class A> auto identity(A&& array) -> decltype(auto) {return std::forward<A>(array);}

template<class A> auto T(A&& arr) -> decltype(auto) {return std::forward<A>(arr).rotated();}  // NOLINT(readability-identifier-naming) : conventional one-letter operation BLAS
template<class A> auto N(A&& arr) -> decltype(auto) {return std::forward<A>(arr)          ;}  // NOLINT(readability-identifier-naming) : conventional one-letter operation BLAS

template<class A, typename D=std::decay_t<A>, typename E=typename D::element_type>
auto conjugated_transposed(A&& arr) -> decltype(auto) {
	return blas::T(blas::conj(std::forward<A>(arr)));
}

// template<class A, typename D=std::decay_t<A>, typename E=typename D::element_type>
// auto conjugated(A&& array) -> decltype(auto) {
//  return blas::conj(std::forward<A>(array));
// }

template<class A>
auto hermitized(A&& array, std::true_type /*true */) -> decltype(auto) {
	return conjugated_transposed(std::forward<A>(array));
}

template<class A>
auto hermitized(A&& array, std::false_type /*false*/) -> decltype(auto) {
	return transposed(std::forward<A>(array));
}

template<class A>
auto hermitized(A&& array) -> decltype(auto) {return conjugated_transposed(std::forward<A>(array));}

namespace operators {

struct H_t {  // NOLINT(readability-identifier-naming) blas naming

	template<class A, std::enable_if_t<std::decay_t<A>::rank::value == 2, int> =0>  // NOLINT(modernize-use-constraints) TODO(correaa) for C++20
	[[nodiscard]] auto operator()(A&& array) const -> decltype(auto) { return hermitized(std::forward<A>(array)); }

	template<class A, 
		std::enable_if_t<std::decay_t<A>::rank::value == 1, int> =0>  // NOLINT(modernize-use-constraints) TODO(correaa) for C++20
	[[deprecated("use blas::C instead of blas::H for conjugated vectors to avoid confusions")]]
	[[nodiscard]] auto operator()(A&& array) const -> decltype(auto) { return blas::conj(std::forward<A>(array)); }

};

inline constexpr H_t H;  // NOLINT(readability-identifier-length) conventional name in BLAS

template<class A, class Op>
auto operator^(A&& array, Op op)
->decltype(op(std::forward<A>(array))) {
	return op(std::forward<A>(array)); }

} // end namespace operators

using operators::H;

template<class A,
	std::enable_if_t<std::decay_t<A>::rank::value == 1, int> =0>  // NOLINT(modernize-use-constraints) TODO(correaa) for C++20
auto C(A&& array) -> decltype(auto) { return blas::conj(std::forward<A>(array)); }  // NOLINT(readability-identifier-naming,readability-identifier-length) : conventional one-letter operation BLAS

template<class A,
	std::enable_if_t<std::decay_t<A>::rank::value == 2, int> =0>  // NOLINT(modernize-use-constraints) TODO(correaa) for C++20
[[deprecated("use blas::H instead of blas::C for conjugated transposed matrices to avoid confusion, use blas::J for only-conjugation of matrices")]]
auto C(A&& array) -> decltype(auto) { return hermitized(std::forward<A>(array)); }  // NOLINT(readability-identifier-naming,readability-identifier-length) : conventional one-letter operation BLAS

template<class A,
	std::enable_if_t<std::decay_t<A>::rank::value == 2, int> =0>  // NOLINT(modernize-use-constraints) for C++20
auto J(A&& array) -> decltype(auto) { return blas::conj(std::forward<A>(array)); }  // NOLINT(readability-identifier-naming,readability-identifier-length) : conventional one-letter operation BLAS

namespace operators {

	template<class A>
	auto operator*(A&& array)
	->decltype(blas::conj(std::forward<A>(array))) {
		return blas::conj(std::forward<A>(array)); }

	template<class A>
	auto operator~(A&& array)
	->decltype(blas::T(std::forward<A>(array))) {
		return blas::T(std::forward<A>(array)); }

} // end namespace operators

} // end namespace boost::multi::blas

#endif
