//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_CONVERT2REAL_H
#define QMCPLUSPLUS_CONVERT2REAL_H

#include <complex>
#include "complex_help.hpp"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "OhmmsPETE/Tensor.h"
#include "OhmmsPETE/OhmmsVector.h"
#include "OhmmsPETE/TinyVector.h"

namespace qmcplusplus
{
/** generic conversion from type T1 to type T2 using implicit conversion
*/
template<typename T1, typename T2, IsReal<T2> = true>
inline void convertToReal(const T1& in, T2& out)
{
  out = static_cast<T2>(in);
}

/** specialization of conversion from complex to real
*/
template<typename T1, typename T2, IsReal<T2> = true>
inline void convertToReal(const std::complex<T1>& in, T2& out)
{
  out = in.real();
}

/* specialization of D-dim vectors
 *
 */
template<typename T1, typename T2, unsigned D>
inline void convertToReal(const TinyVector<T1, D>& in, TinyVector<T2, D>& out)
{
  for (int i = 0; i < D; ++i)
    convertToReal(in[i], out[i]);
}

/** specialization for D tensory*/
template<typename T1, typename T2, unsigned D>
inline void convertToReal(const Tensor<T1, D>& in, Tensor<T2, D>& out)
{
  for (int i = 0; i < D * D; ++i)
    convertToReal(in[i], out[i]);
}

/** generic function to convert arrays
 * @param in starting address of type T1
 * @param out starting address of type T2
 * @param n size of in/out
 */
template<typename T1, typename T2>
inline void convertToReal(const T1* restrict in, T2* restrict out, std::size_t n)
{
  for (int i = 0; i < n; ++i)
    convertToReal(in[i], out[i]);
}

/** specialization for a vector */
template<typename T1, typename T2>
inline void convertToReal(const Vector<T1>& in, Vector<T2>& out)
{
  convertToReal(in.data(), out.data(), in.size());
}

/** specialization for a vector */
template<typename T1, typename T2>
inline void convertToReal(const Matrix<T1>& in, Matrix<T2>& out)
{
  convertToReal(in.data(), out.data(), in.size());
}

} // namespace qmcplusplus
#endif
