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
#ifndef QMCPLUSPLUS_MATRIX_ACCESSOR_HPP
#define QMCPLUSPLUS_MATRIX_ACCESSOR_HPP

namespace qmcplusplus
{
namespace testing
{
enum class MA_ori
{
  ROW,
  COLUMN
};

/** Read only access, for testing!
 *
 *  Makes checkMatrix dependency free
 */
template<typename T>
class MatrixAccessor
{
public:
  using value_type = T;

  MatrixAccessor(const T* data, size_t m, size_t n, MA_ori ori = MA_ori::ROW)
      : m_(m),
        n_(n),
        data_(data),
        ori_(ori){

        };

  T operator()(size_t i, size_t j) const
  {
    if (ori_ == MA_ori::ROW)
      return data_[i * n_ + j];
    else
      return data_[i + j * m_];
  }

  size_t rows() { return m_; }
  size_t cols() { return n_; }

private:
  size_t m_;
  size_t n_;
  const T* data_;
  MA_ori ori_;
};

} // namespace testing
} // namespace qmcplusplus
#endif
