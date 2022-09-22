// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2019-2022 Alfredo A. Correa

#ifndef MULTI_ADAPTORS_BLAS_OPERATIONS_HPP
#define MULTI_ADAPTORS_BLAS_OPERATIONS_HPP

#include "../blas/numeric.hpp"

namespace boost::multi::blas {

template<class M> auto transposed(M const& array) -> decltype(auto){return rotated(array);}

template<class A, typename D=std::decay_t<A>, typename E=typename D::element_type>
auto conjugated_transposed(A&& array) -> decltype(auto) {
	return transposed(blas::conj(std::forward<A>(array)));
}

template<class A> auto identity(A&& array) -> decltype(auto) {return std::forward<A>(array);}

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

template<class A>
auto transposed(A&& array) -> decltype(auto) {return rotated(std::forward<A>(array));}

namespace operators {

[[maybe_unused]] constexpr static struct {

	template<class A, std::enable_if_t<std::decay_t<A>::rank_v == 2, int> =0>
	auto operator()(A&& array) const -> decltype(auto){return hermitized(std::forward<A>(array));}

	template<class A, std::enable_if_t<std::decay_t<A>::rank_v == 1, int> =0>
	[[deprecated("use blas::C instead of blas::H for conjugated vectors to avoid confusions")]]
	auto operator()(A&& array) const -> decltype(auto){return blas::conj(std::forward<A>(array));}

} H;  // NOLINT(readability-identifier-length) conventional name in BLAS

template<class A, class Op>
auto operator^(A&& array, Op op)
->decltype(op(std::forward<A>(array))){
	return op(std::forward<A>(array));}

} // end namespace operators

using operators::H;

template<class A, std::enable_if_t<std::decay_t<A>::rank_v == 1, int> =0>
auto C(A&& array) -> decltype(auto) {return blas::conj(std::forward<A>(array));}  // NOLINT(readability-identifier-naming,readability-identifier-length) : conventional one-letter operation BLAS

template<class A, std::enable_if_t<std::decay_t<A>::rank_v == 2, int> =0>
auto C(A&& array) -> decltype(auto) {return hermitized(std::forward<A>(array));}  // NOLINT(readability-identifier-naming,readability-identifier-length) : conventional one-letter operation BLAS

namespace operators{

	template<class A>
	auto operator*(A&& array)
	->decltype(blas::conj(std::forward<A>(array))){
		return blas::conj(std::forward<A>(array));}

} // end namespace operators

template<class A> auto T(A&& array) -> decltype(auto) {return transposed(std::forward<A>(array));}  // NOLINT(readability-identifier-naming) : conventional one-letter operation BLAS
template<class A> auto N(A&& array) -> decltype(auto) {return identity  (std::forward<A>(array));}  // NOLINT(readability-identifier-naming) : conventional one-letter operation BLAS

} // end namespace boost::multi::blas

#endif
