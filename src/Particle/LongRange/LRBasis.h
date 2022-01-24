//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_LRBASIS_H
#define QMCPLUSPLUS_LRBASIS_H

#include "Particle/ParticleSet.h"
#include "Pools/PooledData.h"
#include "coulomb_types.h"

namespace qmcplusplus
{
/** @ingroup longrange
 *\brief Base-class for long-range breakups.
 *
 * Contains 3 important functions: c(n,k), h(n,r), hintr2(r)
 *  which evaluate the n'th basis function at k or r.
 *  \f[ hintr2_n = \int dr h_n(r) r^2, \f]
 */

struct LRBasis
{
  DECLARE_COULOMB_TYPES

  ///size of the basis elements
  int BasisSize;
  ///Real-space cutoff for short-range part
  mRealType m_rc;

  ///Typedef for the lattice-type. We don't need the full particle-set.
  using ParticleLayout = ParticleSet::ParticleLayout;

  //A copy of the lattice, so that we have the cell-vectors + volume.
  ParticleLayout Lattice;

public:
  //Constructor: copy lattice and init rc
  LRBasis(const ParticleLayout& ref) : m_rc(0.0), Lattice(ref)
  { /*Do nothing*/
  }

  virtual ~LRBasis() = default;

  inline int NumBasisElem() const { return BasisSize; }

  //Real-space basis function + integral: override these
  virtual mRealType h(int n, mRealType r) const = 0;
  virtual mRealType hintr2(int n) const         = 0;
  virtual mRealType dh_dr(int n, mRealType r) const { return 0.0; };
  //k-space basis function: override this
  virtual mRealType c(int m, mRealType k) const = 0;
  //k-space basis function k space derivative.
  virtual mRealType dc_dk(int m, mRealType k) const { return 0.0; };

  //
  // df(m,r) is included for legacy reasons.  Please use dh_dr
  //

  inline mRealType df(int m, mRealType r) const { return dh_dr(m, r); };

  ///$f(r,{tn})$ returns the value of $\sum_n t_n*h_{\alpha n}(r)$
  ///   r is radial position (scalar)
  ///   std::vector<RealType> coefs are the {tn} optimized breakup coefficients.

  /*
 * 
 * name: f: function interpolated by optimized-breakup.
 * @param r The radial coordinate.
 * @param coefs A vector of optimized breakup coefficents.
 * @return The value of $f(r) = \sum_n t_n h_{n \alpha}(r)$
 * 
 */

  inline mRealType f(mRealType r, const std::vector<mRealType>& coefs) const
  {
    mRealType f = 0.0;
    //RealType df = myFunc.df(r, rinv);
    for (int n = 0; n < coefs.size(); n++)
      f += coefs[n] * h(n, r);
    return f;
  }

  /*
 * 
 * name: df_dr
 * @param r The radial coordinate.
 * @param coefs A vector of optimized breakup coefficents.
 * @return $\frac{df}{dr}$.
 * 
 */

  inline mRealType df_dr(mRealType r, const std::vector<mRealType>& coefs) const
  {
    mRealType df = 0.0;
    //RealType df = myFunc.df(r, rinv);
    for (int n = 0; n < coefs.size(); n++)
      df += coefs[n] * dh_dr(n, r);
    return df;
  }

  /*
 * 
 * name: fk
 * @param k |k|
 * @param coefs Optimized breakup coefficients
 * @return The fourier transform of $f(r)$
 * 
 */


  inline mRealType fk(mRealType k, const std::vector<mRealType> coefs) const
  {
    mRealType fk = 0.0;
    for (int n = 0; n < coefs.size(); n++)
      fk += coefs[n] * c(n, k);
    return fk;
  }

  /*
 * 
 * name: dfk_dk
 * @param k |k|
 * @param coefs Optimized breakup coefficients
 * @return $\frac{df_k}{dk}$
 * 
 */

  inline mRealType dfk_dk(mRealType k, const std::vector<mRealType> coefs) const
  {
    mRealType dfk = 0.0;
    for (int n = 0; n < coefs.size(); n++)
      dfk += coefs[n] * dc_dk(n, k);
    return dfk;
  }

  //May need extra functionality when resetting rc. Override this.
  virtual void set_rc(mRealType rc) = 0;
  inline mRealType get_rc() const { return m_rc; }
  inline mRealType get_CellVolume() const { return Lattice.Volume; }
  inline void set_Lattice(const ParticleLayout& ref) { Lattice = ref; }
  inline const ParticleLayout& get_Lattice() const { return Lattice; }
};

} // namespace qmcplusplus

#endif
