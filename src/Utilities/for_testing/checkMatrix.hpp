//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_CHECKMATRIX_HPP
#define QMCPLUSPLUS_CHECKMATRIX_HPP

#include "catch.hpp"

#include <string>
#include <complex>
#include <type_traits>
#include "type_traits/complex_help.hpp"
namespace qmcplusplus
{

template<typename T, IsComplex<T> = true>
bool approxEquality(T val_a, T val_b)
{
  return val_a == ComplexApprox(val_b);
}

template<typename T, IsReal<T> = true>
bool approxEquality(T val_a, T val_b)
{
  return val_a == Approx(val_b);
}

struct CheckMatrixResult
{
  bool result;
  std::string result_message;
};

/** This function checks equality a_mat and b_mat elements
 *  M1, M2 need to have their element type declared M1::value_type
 *         and have an operator(i,j) accessor.
 *  I leave the c++14 template meta programming to insure 
 *  this as an exercise for the reader. Or just enjoy the compiler error.
 *
 *  \param[in] a_mat     - reference matrix, if padded must be identical to b_mat,
 *                         can be a smaller than b_mat in which case it is compared to upper
 *                         left block of b_mat.
 *  \param[in] b_mat     - the matrix to check
 *  \param[in] check_all - if true continue to check matrix elements after failure
 */
template<class M1, class M2>
CheckMatrixResult checkMatrix(M1& a_mat, M2& b_mat, const bool check_all = false)
{
  // This allows use to check a padded b matrix with a nonpadded a
  if (a_mat.rows() > b_mat.rows() || a_mat.cols() > b_mat.cols())
    return {false, "b_mat is too small for a_mat to be a checkable block"};
  std::stringstream error_msg;
  auto matrixElementError = [&error_msg](int i, int j, auto& a_mat, auto& b_mat) {
    error_msg << "checkMatrix found bad element at " << i << ":" << j << "  " << a_mat(i, j) << " != " << b_mat(i, j)
              << '\n';
  };
  bool all_elements_match = true;
  for (int i = 0; i < a_mat.rows(); i++)
    for (int j = 0; j < a_mat.cols(); j++)
    {
      bool approx_equality = approxEquality<typename M1::value_type>(a_mat(i, j), b_mat(i, j));
      if (!approx_equality)
      {
        matrixElementError(i, j, a_mat, b_mat);
        all_elements_match = false;
        if (!check_all)
          return {false, error_msg.str()};
      }
    }
  return {all_elements_match, error_msg.str()};
}

extern template bool approxEquality<float>(float val_a, float val_b);
extern template bool approxEquality<std::complex<float>>(std::complex<float> val_a, std::complex<float> val_b);
extern template bool approxEquality<double>(double val_a, double val_b);
extern template bool approxEquality<std::complex<double>>(std::complex<double> val_a, std::complex<double> val_b);
} // namespace qmcplusplus
#endif
