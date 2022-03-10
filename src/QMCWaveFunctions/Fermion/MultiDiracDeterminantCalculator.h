//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//                    ChangMo Yang, nichthierwohne@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: ChangMo Yang, nichthierwohne@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


/**@file
 */
#ifndef QMCPLUSPLUS_MULTIDIRACDETERMINANTCALCULATOR_H
#define QMCPLUSPLUS_MULTIDIRACDETERMINANTCALCULATOR_H

#include "OhmmsPETE/OhmmsMatrix.h"
#include "Numerics/DeterminantOperators.h"

namespace qmcplusplus
{
/** Function class calculates multi determinant ratio from matrix elements
 *  needs to be the size of your result matrix or larger
 *  includes manual expansions for smaller evaluations
 */
template<typename T>
struct MultiDiracDeterminantCalculator
{
  std::vector<T> M;
  std::vector<int> Pivot;

  void resize(int n)
  {
    M.resize(n * n);
    Pivot.resize(n);
  }

  inline T evaluate(T a11, T a12, T a21, T a22) { return a11 * a22 - a21 * a12; }

  inline T evaluate(T a11, T a12, T a13, T a21, T a22, T a23, T a31, T a32, T a33)
  {
    return (a11 * (a22 * a33 - a32 * a23) - a21 * (a12 * a33 - a32 * a13) + a31 * (a12 * a23 - a22 * a13));
  }

  inline T evaluate(T a11,
                    T a12,
                    T a13,
                    T a14,
                    T a21,
                    T a22,
                    T a23,
                    T a24,
                    T a31,
                    T a32,
                    T a33,
                    T a34,
                    T a41,
                    T a42,
                    T a43,
                    T a44)
  {
    return (a11 * (a22 * (a33 * a44 - a43 * a34) - a32 * (a23 * a44 - a43 * a24) + a42 * (a23 * a34 - a33 * a24)) -
            a21 * (a12 * (a33 * a44 - a43 * a34) - a32 * (a13 * a44 - a43 * a14) + a42 * (a13 * a34 - a33 * a14)) +
            a31 * (a12 * (a23 * a44 - a43 * a24) - a22 * (a13 * a44 - a43 * a14) + a42 * (a13 * a24 - a23 * a14)) -
            a41 * (a12 * (a23 * a34 - a33 * a24) - a22 * (a13 * a34 - a33 * a14) + a32 * (a13 * a24 - a23 * a14)));
  }

  inline T evaluate(T a11,
                    T a12,
                    T a13,
                    T a14,
                    T a15,
                    T a21,
                    T a22,
                    T a23,
                    T a24,
                    T a25,
                    T a31,
                    T a32,
                    T a33,
                    T a34,
                    T a35,
                    T a41,
                    T a42,
                    T a43,
                    T a44,
                    T a45,
                    T a51,
                    T a52,
                    T a53,
                    T a54,
                    T a55)
  {
    return (
        a11 *
            (a22 * (a33 * (a44 * a55 - a54 * a45) - a43 * (a34 * a55 - a54 * a35) + a53 * (a34 * a45 - a44 * a35)) -
             a32 * (a23 * (a44 * a55 - a54 * a45) - a43 * (a24 * a55 - a54 * a25) + a53 * (a24 * a45 - a44 * a25)) +
             a42 * (a23 * (a34 * a55 - a54 * a35) - a33 * (a24 * a55 - a54 * a25) + a53 * (a24 * a35 - a34 * a25)) -
             a52 * (a23 * (a34 * a45 - a44 * a35) - a33 * (a24 * a45 - a44 * a25) + a43 * (a24 * a35 - a34 * a25))) -
        a21 *
            (a12 * (a33 * (a44 * a55 - a54 * a45) - a43 * (a34 * a55 - a54 * a35) + a53 * (a34 * a45 - a44 * a35)) -
             a32 * (a13 * (a44 * a55 - a54 * a45) - a43 * (a14 * a55 - a54 * a15) + a53 * (a14 * a45 - a44 * a15)) +
             a42 * (a13 * (a34 * a55 - a54 * a35) - a33 * (a14 * a55 - a54 * a15) + a53 * (a14 * a35 - a34 * a15)) -
             a52 * (a13 * (a34 * a45 - a44 * a35) - a33 * (a14 * a45 - a44 * a15) + a43 * (a14 * a35 - a34 * a15))) +
        a31 *
            (a12 * (a23 * (a44 * a55 - a54 * a45) - a43 * (a24 * a55 - a54 * a25) + a53 * (a24 * a45 - a44 * a25)) -
             a22 * (a13 * (a44 * a55 - a54 * a45) - a43 * (a14 * a55 - a54 * a15) + a53 * (a14 * a45 - a44 * a15)) +
             a42 * (a13 * (a24 * a55 - a54 * a25) - a23 * (a14 * a55 - a54 * a15) + a53 * (a14 * a25 - a24 * a15)) -
             a52 * (a13 * (a24 * a45 - a44 * a25) - a23 * (a14 * a45 - a44 * a15) + a43 * (a14 * a25 - a24 * a15))) -
        a41 *
            (a12 * (a23 * (a34 * a55 - a54 * a35) - a33 * (a24 * a55 - a54 * a25) + a53 * (a24 * a35 - a34 * a25)) -
             a22 * (a13 * (a34 * a55 - a54 * a35) - a33 * (a14 * a55 - a54 * a15) + a53 * (a14 * a35 - a34 * a15)) +
             a32 * (a13 * (a24 * a55 - a54 * a25) - a23 * (a14 * a55 - a54 * a15) + a53 * (a14 * a25 - a24 * a15)) -
             a52 * (a13 * (a24 * a35 - a34 * a25) - a23 * (a14 * a35 - a34 * a15) + a33 * (a14 * a25 - a24 * a15))) +
        a51 *
            (a12 * (a23 * (a34 * a45 - a44 * a35) - a33 * (a24 * a45 - a44 * a25) + a43 * (a24 * a35 - a34 * a25)) -
             a22 * (a13 * (a34 * a45 - a44 * a35) - a33 * (a14 * a45 - a44 * a15) + a43 * (a14 * a35 - a34 * a15)) +
             a32 * (a13 * (a24 * a45 - a44 * a25) - a23 * (a14 * a45 - a44 * a15) + a43 * (a14 * a25 - a24 * a15)) -
             a42 * (a13 * (a24 * a35 - a34 * a25) - a23 * (a14 * a35 - a34 * a15) + a33 * (a14 * a25 - a24 * a15))));
  }

  /** default implementation of MultiDiracDeterminant::CalculateRatioFromMatrixElements
   *  If you don't have a perfect square below 25 of dots this is what you hit
   *  dots must be a 2n by 2n Matrix
   *  iter must be size n^2
   *  this object must have been resized to n
   */
  template<typename DT>
  using OffloadMatrix = Matrix<DT, OffloadPinnedAllocator<DT>>;
  template<typename ITER>
  inline T evaluate(OffloadMatrix<T>& dots, ITER it, int n)
  {
    typename std::vector<T>::iterator d = M.begin();
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
      {
        //performance through proper iterator indistiquishable from data pointer
        *(d++) = dots(*(it + i), *(it + n + j));
      }
    return Determinant(M.data(), n, n, Pivot.data());
  }
};

} // namespace qmcplusplus
#endif
