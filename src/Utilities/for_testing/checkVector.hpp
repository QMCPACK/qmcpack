//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2025 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//////////////////////////////////////////////////////////////////////////////////////

/** \file
 *  testing utility function to check vector indexed its return value structure
 *  this not intended to be used in actual application code
 */

#ifndef QMCPLUSPLUS_CHECKVECTOR_HPP
#define QMCPLUSPLUS_CHECKVECTOR_HPP

#include "catch.hpp"

#include <string>
#include <complex>
#include <type_traits>
#include <optional>
#include "OhmmsPETE/TinyVector.h"
#include "type_traits/complex_help.hpp"
#include "ApproximateEquality.hpp"
namespace qmcplusplus
{

/** return structure from vector check
 *  For clean use with catch2 CHECKED_ELSE macro
 *  pass/fail --> true/false in result.
 *  result_message should/must contain information in case of fail.
 *  We reserve the right to report information in result message in the case of pass.
 *
 *  DO NOT use result_message presense or contents to determine pass fail.
 *  If this bothers you see
 *  https://isocpp.github.io/CppCoreGuidelines/CppCoreGuidelines#f21-to-return-multiple-out-values-prefer-returning-a-struct
 */
struct CheckVectorResult
{
  /** If matrix check is successful result = true
   *  if one or more elements fails result = false
   */
  bool result;
  /**  always a valid std::string object
   *   result = true --> default constructured string at this time
   *   result = false --> result_message element failure information,
   *   the extent of which is determined by check_all flag to checkVector
   */
  std::string result_message;
};

/** This function checks equality a_vec and b_vec elements
 *  M1, M2 need to have their element type declared M1::value_type
 *         and have an operator(i,j) accessor.
 *
 *  \param[in] a_vec     - reference vector, if padded must be identical to b_vec,
 *                         can be a smaller than b_vec in which case it is compared to upper
 *                         left block of b_vec.
 *  \param[in] b_vec     - the vectorto check
 *  \param[in] check_all - if true continue to check vector elements after failure
 *  \param[in] eps       - add a tolerance for Catch Approx checks. Default to same as in Approx.
 *  The semantics of the return value are discussed above.
 */
template<unsigned D>
CheckVectorResult checkVector(const std::vector<TinyVector<int, D>>& a_vec,
                              const std::vector<TinyVector<int, D>>& b_vec,
                              const bool check_all = false)
{
  // This allows use to check a padded b matrix with a nonpadded a
  if (a_vec.size() > b_vec.size())
    return {false, "b_vec is too small for a_vec to be a checkable segment"};
  std::stringstream result_msg;
  auto vectorElementError = [&result_msg](int i, auto& a_vec, auto& b_vec) {
    result_msg << "checkVector found bad element at " << i << "  " << a_vec[i] << " != " << b_vec[i] << '\n';
  };
  bool all_elements_match = true;
  for (int i = 0; i < a_vec.size(); i++)
  {
    bool equality = a_vec[i] == b_vec[i];
    if (!equality)
    {
      vectorElementError(i, a_vec, b_vec);
      all_elements_match = false;
      if (!check_all)
        return {false, result_msg.str()};
    }
  }
  return {all_elements_match, result_msg.str()};
}

template<class M1, class M2>
CheckVectorResult checkVector(const M1& a_vec,
                              const M2& b_vec,
                              const bool check_all            = false,
                              std::optional<const double> eps = std::nullopt)
{
  // This allows use to check a padded b matrix with a nonpadded a
  if (a_vec.size() > b_vec.size())
    return {false, "b_vec is too small for a_vec to be a checkable segment"};
  std::stringstream result_msg;
  auto vectorElementError = [&result_msg](int i, auto& a_vec, auto& b_vec) {
    result_msg << "checkVector found bad element at " << i << "  " << a_vec[i] << " != " << b_vec[i] << '\n';
  };
  bool all_elements_match = true;
  for (int i = 0; i < a_vec.size(); i++)
  {
    bool approx_equality = approxEquality<typename M1::value_type>(a_vec[i], b_vec[i], eps);
    if (!approx_equality)
    {
      vectorElementError(i, a_vec, b_vec);
      all_elements_match = false;
      if (!check_all)
        return {false, result_msg.str()};
    }
  }
  return {all_elements_match, result_msg.str()};
}
} // namespace qmcplusplus

#endif
