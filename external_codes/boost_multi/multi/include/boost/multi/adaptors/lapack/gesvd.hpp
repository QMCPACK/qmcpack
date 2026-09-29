// Copyright 2025 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_MULTI_ADAPTORS_LAPACK_GESVD_HPP
#define BOOST_MULTI_ADAPTORS_LAPACK_GESVD_HPP

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <stdexcept>
#include <string>

#define MLP_INT int32_t  // int64_t  // #define INT int64_t

#include "boost/multi/array.hpp"

extern "C" {
void dgesvd_(char const& jobu, char const& jobvt, int const& mm, int const& nn, double* aa, int const& lda, double* ss, double* uu, int const& ldu, double* vt, int const& ldvt, double* work, int const& lwork, int& info);  // NOLINT_INT(readability-identifier-naming)
}

namespace boost::multi::lapack {

template<class Alloc, class AArray2D, class UArray2D, class SArray1D, class VTArray2D>
void gesvd(AArray2D&& AA, UArray2D& UU, SArray1D& ss, VTArray2D& VV, Alloc alloc) {
	assert(AA.size() == UU.size());
	assert((~AA).size() == VV.size());
	assert(ss.size() == std::min(UU.size(), VV.size()));

	assert((~AA).stride() == 1);
	assert(ss.stride() == 1);
	assert((~UU).stride() == 1);
	assert((~VV).stride() == 1);

	MLP_INT info;   // NOLINT(cppcoreguidelines-init-variables) init by function
	double  dwork;  // NOLINT(cppcoreguidelines-init-variables) init by function

	dgesvd_(
		'A' /*all left*/, 'A' /*all right*/,
		static_cast<MLP_INT>(VV.size()), static_cast<MLP_INT>(UU.size()),
		AA.base(), static_cast<MLP_INT>(AA.stride()),
		ss.base(),
		VV.base(), static_cast<MLP_INT>(VV.stride()),
		UU.base(), static_cast<MLP_INT>(UU.stride()),
		&dwork, -1, info
	);
	if(info != 0) {
		throw std::runtime_error("Error in DGESVD work estimation, info: " + std::to_string(info));
	}

	auto const  lwork = static_cast<MLP_INT>(dwork);
	auto* const work  = alloc.allocate(lwork);

	dgesvd_(
		'A' /*all left*/, 'A' /*all right*/,
		static_cast<MLP_INT>(VV.size()), static_cast<MLP_INT>(UU.size()),
		AA.base(), static_cast<MLP_INT>(AA.stride()),
		ss.base(),
		VV.base(), static_cast<MLP_INT>(VV.stride()),
		UU.base(), static_cast<MLP_INT>(UU.stride()),
		work, lwork, info
	);
	alloc.deallocate(work, lwork);

	if(info != 0) {
		throw std::runtime_error("Error in DGESVD computation, info: " + std::to_string(info));
	}

	(void)std::forward<AArray2D>(AA);
}

template<
	template<typename> class AllocT = std::allocator,
	class AArray2D, class UArray2D, class SArray1D, class VTArray2D,
	class Alloc = AllocT<typename std::decay_t<AArray2D>::element>>
void gesvd(AArray2D&& AA, UArray2D&& UU, SArray1D&& ss, VTArray2D&& VV) {
	return gesvd(std::forward<AArray2D>(AA), std::forward<UArray2D>(UU), std::forward<SArray1D>(ss), std::forward<VTArray2D>(VV), Alloc{});
}

template<class Array2D, typename ElementType = typename Array2D::element>
auto gesvd(Array2D const& AA) {
	auto AA_copy = AA;
	auto ret     = std::make_tuple(
        ::boost::multi::array<ElementType, 2>({AA.size(), AA.size()}),             // UU  // Right singular vectors
        ::boost::multi::array<ElementType, 1>(std::min(AA.size(), (~AA).size())),  // ss  // Singular values
        ::boost::multi::array<ElementType, 2>({(~AA).size(), (~AA).size()})        // VV  // Left singular vectors
    );
	gesvd(AA_copy, std::get<0>(ret), std::get<1>(ret), std::get<2>(ret));
#if defined(_MULTI_USING_LAPACK_MKL)
	std::get<2>(ret) = +~std::get<2>(ret);  // MKL returns the transpose of the right singular vectors
#endif
	return ret;
}

}  // end namespace boost::multi::lapack

#endif  // BOOST_MULTI_ADAPTORS_LAPACK_GESVD_HPP
