//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_BAREKINETICHELPER_H
#define QMCPLUSPLUS_BAREKINETICHELPER_H

#include "ParticleBase/ParticleAttribOps.h"

namespace qmcplusplus
{
/** compute real(laplacian)
 */
template<typename T, unsigned D>
inline T laplacian(const TinyVector<T, D>& g, T l)
{
  return dot(g, g) + l;
}

/** specialization of laplacian with complex g & l
 */
template<typename T, unsigned D>
inline T laplacian(const TinyVector<std::complex<T>, D>& g, const std::complex<T>& l)
{
  return l.real() + OTCDot<T, T, D>::apply(g, g);
}

/** Convenience function to compute \f$\Re( \nabla^2_i \partial \Psi_T/\Psi_T)\f$
 * @param g OHMMS_DIM dimensional vector for \f$\nabla_i \ln \Psi_T\f$ .  
 * @param l A number, representing \f$\nabla^2_i \ln \Psi_T\f$ .
 * @param gg OHMMS_DIM dimensional vector containing \f$\nabla_i \partial \ln \Psi_T\f$ . 
 * @param gl A number, representing \f$\nabla^2_i \partial \ln \Psi_T\f$
 * @param ideriv A number, representing \f$\partial \ln \Psi_T\f$
 *
 * @return A number corresponding to \f$\Re( \nabla^2_i \partial \Psi_T/\Psi_T)\f$
 */

template<typename T, unsigned D>
inline T dlaplacian(const TinyVector<T, D>& g, const T l, const TinyVector<T, D>& gg, const T gl, const T ideriv)
{
  return gl + l * ideriv + 2.0 * dot(g, gg) + dot(g, g) * ideriv;
}

template<typename T, unsigned D>
inline T dlaplacian(const TinyVector<std::complex<T>, D>& g,
                    const std::complex<T> l,
                    const TinyVector<std::complex<T>, D>& gg,
                    const std::complex<T> gl,
                    const std::complex<T> ideriv)
{
  std::complex<T> l_times_ideriv                = l * ideriv;
  TinyVector<std::complex<T>, D> g_times_ideriv = g * ideriv;

  return gl.real() + l_times_ideriv.real() + 2.0 * OTCDot<T, T, D>::apply(g, gg) +
      OTCDot<T, T, D>::apply(g, g_times_ideriv);
}

} // namespace qmcplusplus
#endif
