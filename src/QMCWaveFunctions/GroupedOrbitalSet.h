//////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Jeongnim Kim
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
/* @file GroupedOrbialSet.h
 * @brief Declaration/definition of GroupedOrbitalSet
 */
#ifndef QMCPLUSPLUS_GROUPED_SPORBITALSET_H
#define QMCPLUSPLUS_GROUPED_SPORBITALSET_H
#include "QMCWaveFunctions/SPOSetBase.h"

namespace qmcplusplus
{

/** A derived class from SPOSetBase
 *
 * Handles multiple groups of Single-particle orbitals.
 * @example GroupedOrbitalSet<TricubicBsplineSet<T> >
 */
template<class OT>
struct GroupedOrbitalSet: public SPOSetBase
{

  ///the type of single-particle orbtials
  typedef vector<OT*>             SPOContainer_t;
  typedef typename OT::value_type value_type;

  SPOContainer_t Phi;

  ///constructor
  GroupedOrbitalSet(int norbs=0)
  {
    setOrbitalSetSize(norbs);
  }

  /**add a single-particle orbital */
  int add(OT* afunction)
  {
    Phi.push_back(afunction);
    return Phi.size()-1;
  }

  void setOrbitalSetSize(int norbs)
  {
    if(norbs == OrbitalSetSize )
      return;
    OrbitalSetSize=norbs;
    BasisSetSize=norbs;
  }

  void resetParameters(const opt_variables_type& active)
  {
    for(int i=0; i<Phi.size(); i++)
      Phi[i]->resetParameters(active);
  }

  void resetTargetParticleSet(ParticleSet& P) { }

  void evaluate(const ParticleSet& P, int iat, ValueVector_t& psi)
  {
    for(int m=0; m<Phi.size(); m++)
      Phi[m]->evaluate(P.R[iat],psi);
  }

  void evaluate(const ParticleSet& P, int iat
                , ValueVector_t& psi, GradVector_t& dpsi, ValueVector_t& d2psi)
  {
    for(int m=0; m<Phi.size(); m++)
      Phi[m]->evaluate(P.R[iat],psi,dpsi,d2psi);
  }

  void evaluate(const ParticleSet& P, int first, int last
                , ValueMatrix_t& logdet, GradMatrix_t& dlogdet, ValueMatrix_t& d2logdet)
  {
    for(int iat=first, i=0; iat<last; iat++,i++)
      for(int m=0; m<Phi.size(); m++)
        Phi[m]->evaluate(P.R[iat],i,logdet,dlogdet,d2logdet);
  }

  void evaluate_notranspose(const ParticleSet& P, int first, int last
                            , ValueMatrix_t& logdet, GradMatrix_t& dlogdet, ValueMatrix_t& d2logdet)
  {
    APP_ABORT("Need specialization of evaluate_notranspose() for grad_grad_logdet. \n");
  }


  void evaluate_notranspose(const ParticleSet& P, int first, int last
                            , ValueMatrix_t& logdet, GradMatrix_t& dlogdet, HessMatrix_t& grad_grad_logdet)
  {
    APP_ABORT("Need specialization of GroupedOrbitalSet::evaluate_notranspose() for grad_grad_logdet. \n");
  }


};
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/
