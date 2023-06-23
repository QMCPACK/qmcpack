////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License.  See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
////////////////////////////////////////////////////////////////////////////////

#ifndef AFQMC_ARRAY_SIZE_HELP_HPP
#define AFQMC_ARRAY_SIZE_HELP_HPP

#include <array>

namespace qmcplusplus
{
namespace afqmc
{
template<class T> auto generic_sizes(T const& A)
->decltype(std::array<std::size_t, 2>{A.size(0), A.size(1)}) {
	return std::array<std::size_t, 2>{A.size(0), A.size(1)}; }

template<typename T, boost::multi::dimensionality_type DIM, class Alloc>
auto generic_sizes(boost::multi::static_array<T, DIM, Alloc> const& A)
->decltype(A.sizes()) {
	return A.sizes(); }

template<typename T, boost::multi::dimensionality_type DIM, typename ElementPtr>
auto generic_sizes(boost::multi::basic_array<T, DIM, ElementPtr> const& A)
->decltype(A.sizes()) {
	return A.sizes(); }
} // namespace afqmc

} // namespace qmcplusplus

#endif
