//////////////////////////////////////////////////////////////////
// (c) Copyright 2006-  by Kris Delaney and Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
/** @file LRHandlerTemp.h
 * @brief Define a LRHandler with two template parameters
 */
#ifndef QMCPLUSPLUS_EWALD_HANLDER_H
#define QMCPLUSPLUS_EWALD_HANLDER_H

#include "LongRange/LRHandlerBase.h"

namespace qmcplusplus
{

/* LR breakup for the standard Ewald method
 *
 * Quasi-2D Ewald method : J. Phys.: Condens. Matter 16, 891 (2004)
 * http://iopscience.iop.org/0953-8984/16/6/017/
 * Note that \f$ \simga \rightarrow 1/\sigma\f$
 * It is possible to use 3D Ewald but for the bulk system, the optimal breakup method
 * is used.
 */
class EwaldHandler: public LRHandlerBase
{

public:
  ///type of supercell
  int SuperCellEnum;
  /// Related to the Gaussian width: \f$ v_l = v(r)erf(\sigma r)\f$
  RealType Sigma;
  ///Volume of the supercell
  RealType Volume;
  ///Area of the supercell: always z is the slab direction
  RealType Area;
  /** Define prefactors for the mixed boundary conditions
   *
   * For quasi-2D (see Appendix A of JPC)
   * PreFactors[0] = \f$ \frac{2\pi}{A}\f$
   * PreFactors[1] = \f$ \frac{2\pi}{A}\frac{1}{\sigma\pi}\f$
   * PreFactors[2] = \f$ \frac{2\pi}{A}\frac{1}{\sigma\pi}\f$
   * PreFactors[3] = \f$ 2\frac{\sqrt{\pi}}{A*\sigma}-\frac{2\pi}{k*A} erfc(\frac{k}{2\sigma}\f$
   */
  TinyVector<RealType,4> PreFactors;
  ///store |k|
  vector<RealType> kMag;
  /// Constructor
  EwaldHandler(ParticleSet& ref, RealType kc_in=-1.0)
    : LRHandlerBase(kc_in)
  {
    Sigma=LR_kc=ref.Lattice.LR_kc;
  }

  /** "copy" constructor
   * @param aLR LRHandlerTemp
   * @param ref Particleset
   *
   * Copy the content of aLR
   * References to ParticleSet or ParticleLayoutout_t are not copied.
   */
  EwaldHandler(const EwaldHandler& aLR, ParticleSet& ref);

  LRHandlerBase* makeClone(ParticleSet& ref)
  {
    return new EwaldHandler(*this,ref);
  }

  void initBreakup(ParticleSet& ref);

  void Breakup(ParticleSet& ref, RealType rs_in)
  {
    initBreakup(ref);
  }

  void resetTargetParticleSet(ParticleSet& ref) { }

  inline RealType evaluate(RealType r, RealType rinv)
  {
    return erfc(r*Sigma)*rinv;
  }

  /** evaluate the contribution from the long-range part for for spline
   */
  inline RealType evaluateLR(RealType r)
  {
    return -erf(r*Sigma)/r;
  }

  inline RealType evaluateSR_k0()
  {
    return 0.0;
  }

  inline RealType evaluateLR_r0()
  {
    return 2.0*Sigma/std::sqrt(M_PI)+PreFactors[3];
  }

  /**  evaluate the first derivative of the short range part at r
   *
   * @param r  radius
   * @param rinv 1/r
   */
  inline RealType srDf(RealType r, RealType rinv)
  {
    return 0.0;
  }

  void fillFk(KContainer& KList);

  /** evaluate k-dependent
   */
  RealType evaluate_slab(RealType z, const vector<int>& kshell
                         , const ComplexType* restrict rk1, const ComplexType* restrict rk2);

  /** evaluate k=0 term at z
   * @param z distance in the slab direction
   * @param zp z*Sigma
   * @return \f$X*z*erf(z*Sgima) + Y*exp^{-z^2*Simga^2}]\f$
   *
   * Here \f$ X=\frac{2\pi}{A}\f$ and \f$ Y=\frac{2\sqrt{\pi}}{A*Simga}\f$
   */
  inline RealType SlabFunc0(RealType z, RealType zp)
  {
    return PreFactors[0]*z*erf(zp)+PreFactors[1]*std::exp(-zp*zp);
  }

  /** evaluate k!=0 term at z
   * @param ks index of this kshell
   * @param z distance in the slab direction
   * @param zp z*Sigma
   */
  inline RealType SlabFuncK(int ks, RealType z, RealType zp)
  {
    RealType expkz=std::exp(kMag[ks]*z);
    RealType kz0=PreFactors[2]*kMag[ks];//could save this
    return erfc(kz0-zp)/expkz+erfc(kz0+zp)*expkz;
  }
};
}
#endif
/***************************************************************************
 * $RCSfile$   $Author: bkclark $
 * $Revision: 3798 $   $Date: 2009-04-29 00:38:29 -0500 (Wed, 29 Apr 2009) $
 * $Id: EwaldHandler.h 3798 2009-04-29 05:38:29Z bkclark $
 ***************************************************************************/
