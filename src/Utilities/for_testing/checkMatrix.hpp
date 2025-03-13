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
 *  testing utility function to check matrix indexed data and its return value structurep
 *  this not intended to be used in actual application code
 */

#ifndef QMCPLUSPLUS_CHECKMATRIX_HPP
#define QMCPLUSPLUS_CHECKMATRIX_HPP

#include "catch.hpp"

#include <string>
#include <complex>
#include <type_traits>
#include <optional>
#include "type_traits/complex_help.hpp"
#include "ApproximateEquality.hpp"

namespace qmcplusplus
{

/** return structure from matrix check
 *  For easy use with catch2 CHECKED_ELSE macro
 *  pass/fail --> true/false in result.
 *  result_message should/must contain information in case of fail.
 *  We reserve the right to report information in result message in the case of pass.
 *
 *  DO NOT use result_message presense or contents to determine pass fail.
 *  If this bothers you see
 *  https://isocpp.github.io/CppCoreGuidelines/CppCoreGuidelines#f21-to-return-multiple-out-values-prefer-returning-a-struct
 */
struct CheckMatrixResult
{
  /** If matrix check is successful result = true
   *  if one or more elements fails result = false
   */
  bool result;
  /**  always a valid std::string object
   *   result = true --> default constructured string at this time
   *   result = false --> result_message element failure information, the extent of which is determined by check_all flag to checkMatrix
   */
  std::string result_message;
};

/** This function checks equality a_mat and b_mat elements
 *  M1, M2 need to have their element type declared M1::value_type
 *         and have an operator(i,j) accessor.
 *
 *  \param[in] a_mat     - reference matrix, if padded must be identical to b_mat,
 *                         can be a smaller than b_mat in which case it is compared to upper
 *                         left block of b_mat.
 *  \param[in] b_mat     - the matrix to check
 *  \param[in] check_all - if true continue to check matrix elements after failure
 *  \param[in] eps       - add a tolerance for Catch Approx checks. Default to same as in Approx.
 *  The semantics of the return value are discussed above.
 */
template<class M1, class M2>
CheckMatrixResult checkMatrix(const M1& a_mat,
                              const M2& b_mat,
                              const bool check_all            = false,
                              std::optional<const double> eps = std::nullopt)
{
  // This allows use to check a padded b matrix with a nonpadded a
  if (a_mat.rows() > b_mat.rows() || a_mat.cols() > b_mat.cols())
    return {false, "b_mat is too small for a_mat to be a checkable block"};
  std::stringstream result_msg;
  auto matrixElementError = [&result_msg](int i, int j, auto& a_mat, auto& b_mat) {
    result_msg << "checkMatrix found bad element at " << i << ":" << j << "  " << a_mat(i, j) << " != " << b_mat(i, j)
              << '\n';
  };
  bool all_elements_match = true;
  for (int i = 0; i < a_mat.rows(); i++)
    for (int j = 0; j < a_mat.cols(); j++)
    {
      bool approx_equality = approxEquality<typename M1::value_type>(a_mat(i, j), b_mat(i, j), eps);
      if (!approx_equality)
      {
        matrixElementError(i, j, a_mat, b_mat);
        all_elements_match = false;
        if (!check_all)
          return {false, result_msg.str()};
      }
    }
  return {all_elements_match, result_msg.str()};
}

} // namespace qmcplusplus
#endif
