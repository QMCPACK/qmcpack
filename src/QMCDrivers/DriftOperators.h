//////////////////////////////////////////////////////////////////
// (c) Copyright 2005- by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef QMCPLUSPLUS_QMCDRIFTOPERATORS_H
#define QMCPLUSPLUS_QMCDRIFTOPERATORS_H
#include "ParticleBase/ParticleAttribOps.h"
#include "QMCWaveFunctions/OrbitalTraits.h"
namespace qmcplusplus {

  template<class T, unsigned D>
  inline T getDriftScale(T tau, const ParticleAttrib<TinyVector<T,D> >& ga) {
    T vsq=Dot(ga,ga);
    return (vsq<numeric_limits<T>::epsilon())? tau:((-1.0+std::sqrt(1.0+2.0*tau*vsq))/vsq);
  }

  template<class T, unsigned D>
  inline T getDriftScale(T tau, 
      const ParticleAttrib<TinyVector<complex<T>,D> >& ga) {
    T vsq=Dot(ga,ga);
    return (vsq<numeric_limits<T>::epsilon())? tau:((-1.0+std::sqrt(1.0+2.0*tau*vsq))/vsq);
  }

  /** da = scaled(tau)*ga
   * @param tau time step
   * @param qf real quantum forces
   * @param drift drift
   */
  template<class T, unsigned D>
  inline void setScaledDrift(T tau, 
      const ParticleAttrib<TinyVector<T,D> >& qf,
      ParticleAttrib<TinyVector<T,D> >& drift) {
    T s = getDriftScale(tau,qf);
    PAOps<T,D>::scale(s,qf,drift);
  }

  /** da = scaled(tau)*ga
   * @param tau time step
   * @param qf complex quantum forces
   * @param drift drift
   */
  template<class T, unsigned D>
  inline void setScaledDrift(T tau, 
      const ParticleAttrib<TinyVector<complex<T>,D> >& qf,
      ParticleAttrib<TinyVector<T,D> >& drift) {
    T s = getDriftScale(tau,qf);
    PAOps<T,D>::scale(s,qf,drift);
  }

  template<class T, unsigned D>
  inline void assignDrift(T s, 
      const ParticleAttrib<TinyVector<T,D> >& ga,
      ParticleAttrib<TinyVector<T,D> >& da) {
    PAOps<T,D>::scale(s,ga,da);
  }

  template<class T, unsigned D>
  inline void assignDrift(T s, 
      const ParticleAttrib<TinyVector<complex<T>,D> >& ga,
      ParticleAttrib<TinyVector<T,D> >& da) {
    PAOps<T,D>::scale(s,ga,da);
  }
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
