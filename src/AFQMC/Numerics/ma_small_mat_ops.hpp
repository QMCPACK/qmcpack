//////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License.  See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by:
//
// File created by:
// Miguel A. Morales, moralessilva2@llnl.gov
//    Lawrence Livermore National Laboratory
////////////////////////////////////////////////////////////////////////////////

#ifndef MA_SMALL_MAT_OPS_HPP
#define MA_SMALL_MAT_OPS_HPP

namespace ma
{
template<typename T>
inline static T D2x2(T const a11, T const a12, T const a21, T const a22)
{
  return a11 * a22 - a21 * a12;
}

template<typename T>
inline static T I2x2(T const a11, T const a12, T const a21, T const a22, T* M)
{
  T det = a11 * a22 - a21 * a12;
  M[0]  = a22 / det;
  M[1]  = -a12 / det;
  M[2]  = -a21 / det;
  M[3]  = a11 / det;
  return det;
}

template<typename T, class Mat>
inline static T I2x2(T const a11, T const a12, T const a21, T const a22, Mat& M)
{
  T det   = a11 * a22 - a21 * a12;
  M[0][0] = a22 / det;
  M[0][1] = -a12 / det;
  M[1][0] = -a21 / det;
  M[1][1] = a11 / det;
  return det;
}

template<typename T>
inline static T D3x3(T const a11,
                     T const a12,
                     T const a13,
                     T const a21,
                     T const a22,
                     T const a23,
                     T const a31,
                     T const a32,
                     T const a33)
{
  return (a11 * (a22 * a33 - a32 * a23) - a21 * (a12 * a33 - a32 * a13) + a31 * (a12 * a23 - a22 * a13));
}

template<typename T>
inline static T I3x3(T const a11,
                     T const a12,
                     T const a13,
                     T const a21,
                     T const a22,
                     T const a23,
                     T const a31,
                     T const a32,
                     T const a33,
                     T* M)
{
  T det = (a11 * (a22 * a33 - a32 * a23) - a21 * (a12 * a33 - a32 * a13) + a31 * (a12 * a23 - a22 * a13));
  M[0]  = (a22 * a33 - a32 * a23) / det;
  M[1]  = (a13 * a32 - a12 * a33) / det;
  M[2]  = (a12 * a23 - a13 * a22) / det;
  M[3]  = (a23 * a31 - a21 * a33) / det;
  M[4]  = (a11 * a33 - a13 * a31) / det;
  M[5]  = (a13 * a21 - a11 * a23) / det;
  M[6]  = (a21 * a32 - a22 * a31) / det;
  M[7]  = (a12 * a31 - a11 * a32) / det;
  M[8]  = (a11 * a22 - a12 * a21) / det;
  return det;
}

template<typename T, class Mat>
inline static T I3x3(T const a11,
                     T const a12,
                     T const a13,
                     T const a21,
                     T const a22,
                     T const a23,
                     T const a31,
                     T const a32,
                     T const a33,
                     Mat& M)
{
  T det   = (a11 * (a22 * a33 - a32 * a23) - a21 * (a12 * a33 - a32 * a13) + a31 * (a12 * a23 - a22 * a13));
  M[0][0] = (a22 * a33 - a32 * a23) / det;
  M[0][1] = (a13 * a32 - a12 * a33) / det;
  M[0][2] = (a12 * a23 - a13 * a22) / det;
  M[1][0] = (a23 * a31 - a21 * a33) / det;
  M[1][1] = (a11 * a33 - a13 * a31) / det;
  M[1][2] = (a13 * a21 - a11 * a23) / det;
  M[2][0] = (a21 * a32 - a22 * a31) / det;
  M[2][1] = (a12 * a31 - a11 * a32) / det;
  M[2][2] = (a11 * a22 - a12 * a21) / det;
  return det;
}

template<typename T>
inline static T D4x4(T const a11,
                     T const a12,
                     T const a13,
                     T const a14,
                     T const a21,
                     T const a22,
                     T const a23,
                     T const a24,
                     T const a31,
                     T const a32,
                     T const a33,
                     T const a34,
                     T const a41,
                     T const a42,
                     T const a43,
                     T const a44)
{
  return (a11 * (a22 * (a33 * a44 - a43 * a34) - a32 * (a23 * a44 - a43 * a24) + a42 * (a23 * a34 - a33 * a24)) -
          a21 * (a12 * (a33 * a44 - a43 * a34) - a32 * (a13 * a44 - a43 * a14) + a42 * (a13 * a34 - a33 * a14)) +
          a31 * (a12 * (a23 * a44 - a43 * a24) - a22 * (a13 * a44 - a43 * a14) + a42 * (a13 * a24 - a23 * a14)) -
          a41 * (a12 * (a23 * a34 - a33 * a24) - a22 * (a13 * a34 - a33 * a14) + a32 * (a13 * a24 - a23 * a14)));
}

template<typename T>
inline static T D5x5(T const a11,
                     T const a12,
                     T const a13,
                     T const a14,
                     T const a15,
                     T const a21,
                     T const a22,
                     T const a23,
                     T const a24,
                     T const a25,
                     T const a31,
                     T const a32,
                     T const a33,
                     T const a34,
                     T const a35,
                     T const a41,
                     T const a42,
                     T const a43,
                     T const a44,
                     T const a45,
                     T const a51,
                     T const a52,
                     T const a53,
                     T const a54,
                     T const a55)
{
  return (a11 *
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

} // namespace ma

#endif
