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

namespace qmcplusplus {

  /* LR breakup for the standard Ewald method
   */
  class EwaldHandler: public LRHandlerBase {

  public:
    /// Related to the Gaussian width: \f$ v_l = v(r)erf(Sigma*r)\f$
    RealType Sigma;
    ///Volume of the supercell
    RealType Volume;

    //Constructor
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
    
    /**  evaluate the first derivative of the short range part at r
     *
     * @param r  radius
     * @param rinv 1/r
     */
    inline RealType srDf(RealType r, RealType rinv) {
      //have to retern the derivative
      return 0.0;
    }
    
    /** evaluate the contribution from the long-range part for for spline
     */
    inline RealType evaluateLR(RealType r) {
      return erf(r*Sigma)/r;
    }

    inline RealType evaluateSR_k0() 
    {
      return 0.0;
    }

    inline RealType evaluateLR_r0()
    {
      return 2.0*Sigma/std::sqrt(M_PI);
    }
    
    void fillFk(KContainer& KList);
  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author: bkclark $
 * $Revision: 3798 $   $Date: 2009-04-29 00:38:29 -0500 (Wed, 29 Apr 2009) $
 * $Id: EwaldHandler.h 3798 2009-04-29 05:38:29Z bkclark $
 ***************************************************************************/
