// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2019-2022 Alfredo A. Correa

#ifndef MULTI_ADAPTORS_BLAS_OPERATIONS_HPP
#define MULTI_ADAPTORS_BLAS_OPERATIONS_HPP

#include "../blas/numeric.hpp"

namespace boost::multi::blas {

template<class M> auto transposed(M const& m) -> decltype(auto){return rotated(m);}

template<class A, typename D=std::decay_t<A>, typename E=typename D::element_type>
auto conjugated_transposed(A&& a) -> decltype(auto){
	return transposed(blas::conj(std::forward<A>(a)));
}

template<class A> auto identity(A&& a) -> decltype(auto){return std::forward<A>(a);}

template<class A>
auto hermitized(A&& a, std::true_type /*true */) -> decltype(auto){
	return conjugated_transposed(std::forward<A>(a));
}

template<class A>
auto hermitized(A&& a, std::false_type /*false*/) -> decltype(auto){
	return transposed(std::forward<A>(a));
}

template<class A>
auto hermitized(A&& a) -> decltype(auto){return conjugated_transposed(std::forward<A>(a));}

template<class A>
auto transposed(A&& a) -> decltype(auto){return rotated(std::forward<A>(a));}

namespace operators{

MAYBE_UNUSED constexpr static struct {

	template<class A, std::enable_if_t<std::decay_t<A>::rank_v == 2, int> =0>
	auto operator()(A&& a) const -> decltype(auto){return hermitized(std::forward<A>(a));}

	template<class A, std::enable_if_t<std::decay_t<A>::rank_v == 1, int> =0>
	[[deprecated("use blas::C instead of blas::H for conjugated vectors to avoid confusions")]]
	auto operator()(A&& a) const -> decltype(auto){return blas::conj(std::forward<A>(a));}

} H;

template<class A, class Op>
auto operator^(A&& a, Op op)
->decltype(op(std::forward<A>(a))){
	return op(std::forward<A>(a));}

} // end namespace operators

using operators::H;

template<class A, std::enable_if_t<std::decay_t<A>::rank_v == 1, int> =0>
auto C(A&& a) -> decltype(auto){return blas::conj(std::forward<A>(a));}  // NOLINT(readability-identifier-naming) : conventional one-letter operation BLAS

template<class A, std::enable_if_t<std::decay_t<A>::rank_v == 2, int> =0>
auto C(A&& a) -> decltype(auto){return hermitized(std::forward<A>(a));}  // NOLINT(readability-identifier-naming) : conventional one-letter operation BLAS

namespace operators{

	template<class A>
	auto operator*(A&& a)
	->decltype(blas::conj(std::forward<A>(a))){
		return blas::conj(std::forward<A>(a));}

} // end namespace operators

template<class A> auto T(A&& a) -> decltype(auto){return transposed(std::forward<A>(a));}  // NOLINT(readability-identifier-naming) : conventional one-letter operation BLAS
template<class A> auto N(A&& a) -> decltype(auto){return identity  (std::forward<A>(a));}  // NOLINT(readability-identifier-naming) : conventional one-letter operation BLAS

} // end namespace boost::multi::blas

#endif
