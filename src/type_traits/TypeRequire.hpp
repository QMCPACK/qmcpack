//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2018 QMCPACK developers
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//                    Ye Luo, yeluo@anl.gov, Argonne National Lab
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Lab
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_TYPEREQUIRE_HPP
#define QMCPLUSPLUS_TYPEREQUIRE_HPP
/** @file
 *  @brief RequireTypeTraitWhatevers used in more than one location in the code
 *  
 *  These fall into a gap in the naming conventions and perhaps shouldn't be in the
 *  qmcplusplus namespace
 */

#include <type_traits>
#include <iterator>

namespace qmcplusplus
{
template<typename InputIterator>
using RequireInputIterator =
    typename std::enable_if<std::is_convertible<typename std::iterator_traits<InputIterator>::iterator_category,
                                                std::input_iterator_tag>::value>::type;

// template<typename WrapsUniquePtr>
// using RequireUnderlyingUniquePtr =
//     typename std::enable_if<std::is_same<typename WrapsUniquePtr::value_type, typename std::unique_ptr<typename WrapsUniquePtr::value_type>>>;

/** Check rvo optimization actually occurs, link will fail if it doesn't, Only for testing
 *
 *  See: https://stackoverflow.com/questions/35736568/is-there-a-way-to-check-if-rvo-was-applied
 */
template<typename T>
struct force_rvo : public T
{
  force_rvo() {}
  using T::T;
  force_rvo(const force_rvo&);
  force_rvo(force_rvo&&) noexcept;
};


} // namespace qmcplusplus

#endif
