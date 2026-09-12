// Copyright 2020-2026 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_MULTI_ADAPTORS_LAPACK_SYEV_HPP
#define BOOST_MULTI_ADAPTORS_LAPACK_SYEV_HPP

#include "boost/multi/adaptors/blas/filling.hpp"
#include "boost/multi/adaptors/lapack/core.hpp"
#include "boost/multi/array.hpp"

#include <cassert>

namespace boost {
namespace multi {
namespace lapack {

// no `using blas::filling;` here: it would conflict with lapack::filling (lapack/filling.hpp)
// when this header is included together with lapack/potrf.hpp; blas::filling is used qualified

using ::core::syev;

/// Computes eigenvalues and eigenvectors of a real symmetric 2D array, in place, with a caller-provided workspace.
///
/// On successful return the rows of `a` are overwritten with the orthonormal eigenvectors
/// and `w` holds the corresponding eigenvalues in ascending order.
/// The computation runs on the CPU (LAPACK `?syev`) or on the GPU
/// (cuSOLVER/hipSOLVER `?syevd`, see `lapack/cusolver.hpp`) depending on the element
/// pointer type of the arrays; the call is a compile-time error for unsupported pointer/element combinations.
///
/// @tparam Array2D 2D array or subarray type (elements of `float` or `double`)
/// @tparam Array1D 1D array type for the eigenvalues
/// @tparam Array1DW 1D array type for the workspace
/// @param uplo triangle of `a` that is referenced, in the C++ (row-major) view:
///        `blas::filling::upper` means the data is stored on and above the diagonal
/// @param a square symmetric matrix; contiguous in exactly one of its two dimensions; overwritten with eigenvectors (as rows)
/// @param w destination for the eigenvalues; `w.size() == a.size()`, contiguous (`stride() == 1`)
/// @param work scratch array of at least `max(1, 3*a.size() - 1)` contiguous elements, allocated where `a`'s elements live
/// @return a view of `a` restricted to the leading block that converged,
///         `a({0, n - k}, {0, n - k})` with `k` the LAPACK `INFO` failure count;
///         its `.size()` equals `a.size()` exactly when all eigenvalues converged
/// @note An empty `a` returns an empty view; the eigenvector rows follow the row-major convention
///       (each row of the result is one eigenvector), matching `A = Vᵀ·diag(w)·V`
template<class Array2D, class Array1D, class Array1DW>
auto syev(blas::filling uplo, Array2D&& a, Array1D&& w, Array1DW&& work)
	-> decltype(syev('V', uplo == blas::filling::upper ? 'L' : 'U', a.size(), a.base(), a.stride(), w.base(), work.base(), work.size(), std::declval<int&>()), a({0L, 1L}, {0L, 1L})) {
	assert(work.size() >= std::max<decltype(a.size())>(1, 3 * a.size() - 1));
	assert(a.size() == w.size());
	assert(w.stride() == 1);
	assert(work.stride() == 1);

	if(a.size() == 0)
		return std::forward<Array2D>(a)();

	int info = -1;

	if(a.rotated().stride() == 1) {
		syev('V', uplo == blas::filling::upper ? 'L' : 'U', a.size(), a.base(), a.stride(), w.base(), work.base(), work.size(), info);
	} else if(a.stride() == 1) {
		syev('V', uplo == blas::filling::upper ? 'U' : 'L', a.size(), a.base(), a.rotated().stride(), w.base(), work.base(), work.size(), info);
	} else {
		assert(0);
	}  // case not contemplated by lapack

	if(info < 0) {
		assert(0);
	}  // bad argument

	return std::forward<Array2D>(a)({0, size(a) - info}, {0, size(a) - info});
}

/// Computes eigenvalues and eigenvectors of a real symmetric 2D array, in place, allocating the workspace internally.
///
/// Same as the workspace-taking overload, with a scratch array of `max(1, 3*a.size() - 1)`
/// elements allocated through `w`'s allocator (so for device arrays the workspace is a device array).
///
/// @param uplo triangle of `a` that is referenced, in the C++ (row-major) view
/// @param a square symmetric matrix; overwritten with the eigenvectors (as rows)
/// @param w destination for the eigenvalues, in ascending order
/// @return a view of `a` spanning the leading block that converged (full size on success)
template<class Array2D, class Array1D, class Array1DW = typename std::decay_t<Array1D>::decay_type>
auto syev(blas::filling uplo, Array2D&& a, Array1D&& w)
	-> decltype(syev(uplo, std::forward<Array2D>(a), std::forward<Array1D>(w), Array1DW(std::max<decltype(size(a))>(1, 3 * size(a) - 1), w.get_allocator()))) {
	return syev(uplo, std::forward<Array2D>(a), std::forward<Array1D>(w), Array1DW(std::max<decltype(size(a))>(1, 3 * size(a) - 1), w.get_allocator()));
}  // TODO(correaa) obtain automatic size from lapack info routine

// /// Computes eigenvalues and eigenvectors of a real symmetric 2D array, without modifying it.
// ///
// /// @param uplo triangle of `a` that is referenced, in the C++ (row-major) view
// /// @param a square symmetric matrix; not modified
// /// @param w destination for the eigenvalues, in ascending order
// /// @return a newly allocated array holding the eigenvectors as rows
// /// @note The result must be used (the input is `const`, so discarding it discards the whole computation)
// template<class Array2D, class Array1D>
// [[nodiscard]]  // "because input array is const, output gives eigenvectors"
// auto syev(blas::filling uplo, Array2D const& a, Array1D&& w) -> typename Array2D::decay_type {
// 	auto ret = a.decay();
// 	auto l   = syev(uplo, ret, std::forward<Array1D>(w));
// 	if(size(l) != size(a))
// 		assert(0);  // failed
// 	return ret;
// }

/// Computes eigenvalues and eigenvectors of a real symmetric 2D array in place, returning the eigenvalues.
///
/// @param uplo triangle of `a` that is referenced, in the C++ (row-major) view
/// @param a square symmetric matrix; overwritten with the eigenvectors (as rows)
/// @return a newly allocated 1D array (using `a`'s allocator) with the eigenvalues in ascending order
template<class Array2D>
[[nodiscard]]  // "because input array is const, output gives eigenvalues"
auto syev(blas::filling uplo, Array2D&& a) {
	multi::array<typename std::decay_t<Array2D>::element, 1, decltype(get_allocator(a))> eigenvalues(a.size(), a.get_allocator());
	syev(uplo, std::forward<Array2D>(a), eigenvalues);
	return eigenvalues;
}

/// Computes eigenvalues and eigenvectors of a real symmetric 2D array, without modifying it, returning both.
///
/// Usable with structured bindings: `auto const [V, w] = syev(blas::filling::upper, A);`
///
/// @param uplo triangle of `a` that is referenced, in the C++ (row-major) view
/// @param a square symmetric matrix; not modified
/// @return an aggregate with members `eigenvectors` (2D, eigenvectors as rows) and
///         `eigenvalues` (1D, ascending), both allocated with `a`'s allocator
template<class Array2D>
[[nodiscard]]  // "because input array is const, output gives a structured binding of eigenvectors and eigenvactor"
auto syev(blas::filling uplo, Array2D const& a) {
	struct {
		typename Array2D::decay_type eigenvectors;
		typename Array2D::value_type eigenvalues;
	} ret{a, typename Array2D::value_type(size(a), get_allocator(a))};
	auto&& l = syev(uplo, ret.eigenvectors, ret.eigenvalues);
	assert(size(l) == size(a));
	return ret;
}

}  // namespace lapack
}  // namespace multi
}  // namespace boost
#endif  // BOOST_MULTI_ADAPTORS_LAPACK_SYEV_HPP
