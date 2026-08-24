// Copyright 2019-2024 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_MULTI_ADAPTORS_LAPACK_POTRF_HPP
#define BOOST_MULTI_ADAPTORS_LAPACK_POTRF_HPP
#pragma once

#include "boost/multi/array.hpp"
#include "boost/multi/detail/config/NODISCARD.hpp"

#include "../lapack/core.hpp"
#include "../lapack/filling.hpp"

#include "../blas/numeric.hpp"

#include<cassert>

namespace boost::multi::lapack {

// using blas::filling;

using ::core::potrf;

template<class Iterator>
BOOST_MULTI_NODISCARD("result has information of order of minor through .size() member")
auto potrf(filling uplo, Iterator first, Iterator last)
->decltype(potrf(static_cast<char>(uplo), typename std::iterator_traits<Iterator>::difference_type{}, first.base(), stride(first), std::declval<int&>()), Iterator{})
{
	assert( stride(first) == stride(last) );
	assert( first->stride() == 1 );
//  auto lda = stride(first);

	int info;  // NOLINT(cppcoreguidelines-init-variables)
	potrf(static_cast<char>(uplo), std::distance(first, last), first.base(), stride(first), info);

	assert( info >= 0 );
	// if(info > 0) {std::cerr << "warning minor of order " << info << " is not possitive\n";}
	return info==0?last:first + info - 1;
}

template<class A2D>
BOOST_MULTI_NODISCARD("result has information of order of minor through .size() member")
auto potrf(filling uplo, A2D&& A)  // NOLINT(readability-identifier-length) conventional lapack name
->decltype(potrf(uplo, begin(A), end(A)), A({0, 1}))
{
	using lapack::flip;

	if(stride(A) == 1) {
		auto last = potrf(flip(uplo), A.rotated().begin(), A.rotated().end());
		using std::distance;
		return A({0, distance(A.rotated().begin(), last)}, {0, distance(A.rotated().begin(), last)});
	}

	auto last = potrf(uplo, begin(A), end(A));

	using std::distance;
	return std::forward<A2D>(A)({0, distance(begin(A), last)});  // , {0, distance(begin(A), last-1)});
}

template<class A>
struct hermitic_t : private A {
	using underlying_type = A;

	auto underlying() const & -> underlying_type const& {return *this;}
	auto underlying()       & -> underlying_type      & {return *this;}
	auto underlying()      && -> underlying_type     && {return std::move(*this);}

 private:
	lapack::filling side_;

 public:
	auto side() const {return side_;}

	hermitic_t(A const& a, lapack::filling side) : A{a}, side_{side} {}  // NOLINT(readability-identifier-length) conventional lapack name
	using A::size;
};

template<class A> auto hermitic(lapack::filling side, A&& a)  // NOLINT(readability-identifier-length) conventional lapack name
-> hermitic_t<std::decay_t<decltype(std::declval<A>()())>> {
	return {std::forward<A>(a)(), side};
}

template<class HA>
BOOST_MULTI_NODISCARD("result is returned because third argument is const")
auto potrf(HA&& ha) -> decltype(auto) {
	return hermitic(ha.side, potrf(ha.side, std::forward<HA>(ha).underlying()));  // static_cast<typename HA::underlying_type&>(ha)));
}

// orthonormalize rows
template<class A> auto onrm(A&& a, filling uplo /*= filling::upper*/)  // NOLINT(readability-identifier-length) conventional lapack name
->decltype(trsm(flip(uplo), hermitized(potrf(uplo, herk(uplo, a))), std::forward<A>(a))) { assert(size(a) <= size(rotated(a)));
	return trsm(flip(uplo), hermitized(potrf(uplo, herk(uplo, a))), std::forward<A>(a)); }

template<class A, class B> auto onrm(A&& a, B&& buffer, filling uplo /* = filling::upper*/)  // NOLINT(readability-identifier-length) conventional lapack name
->decltype(trsm(flip(uplo), hermitized(potrf(uplo, herk(uplo, a, std::forward<B>(buffer)))), std::forward<A>(a))) { assert(size(a) <= size(rotated(a)));
	return trsm(flip(uplo), hermitized(potrf(uplo, herk(uplo, a, std::forward<B>(buffer)))), std::forward<A>(a)); }

}  // end namespace boost::multi::lapack
#endif
