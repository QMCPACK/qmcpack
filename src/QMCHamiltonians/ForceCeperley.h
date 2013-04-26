//////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Jeongnim Kim and Kris Delaney
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
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
#ifndef QMCPLUSPLUS_FORCE_CEPERLEY_HAMILTONIAN_H
#define QMCPLUSPLUS_FORCE_CEPERLEY_HAMILTONIAN_H
#include "QMCHamiltonians/ForceBase.h"
#include "LongRange/LRCoulombSingleton.h"
#include "Numerics/OneDimGridBase.h"
#include "Numerics/OneDimGridFunctor.h"
#include "Numerics/OneDimCubicSpline.h"

namespace qmcplusplus
{

struct ForceCeperley: public QMCHamiltonianBase, public ForceBase
{

  double Rcut; // parameter: radial distance within which estimator is used
  int m_exp; // parameter: exponent in polynomial fit
  int N_basis; // parameter: size of polynomial basis set
  Matrix<RealType> Sinv; // terms in fitting polynomial
  Vector<double> h; // terms in fitting polynomial
  Vector<double> c; // polynomial coefficients
  // container for short-range force estimator

  ParticleSet::ParticlePos_t forces_ShortRange;

  ForceCeperley(ParticleSet& ions, ParticleSet& elns);

  Return_t evaluate(ParticleSet& P);

  void InitMatrix();

  inline Return_t evaluate(ParticleSet& P, vector<NonLocalData>& Txy)
  {
    return evaluate(P);
  }

  void registerObservables(vector<observable_helper*>& h5list, hid_t gid) const
  {
    registerObservablesF(h5list,gid);
  }

  void addObservables(PropertySetType& plist, BufferType& collectables)
  {
    addObservablesF(plist);
  }

  void setObservables(PropertySetType& plist)
  {
    setObservablesF(plist);
  }

  void resetTargetParticleSet(ParticleSet& P)
  {
  }

  void setParticlePropertyList(PropertySetType& plist, int offset)
  {
    setParticleSetF(plist, offset);
  }
  QMCHamiltonianBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi);

  bool put(xmlNodePtr cur) ;

  bool get(std::ostream& os) const
  {
    os << "Ceperley Force Estimator Hamiltonian: " << pairName;
    return true;
  }

};

}
#endif

/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 2945 $   $Date: 2008-08-05 10:21:33 -0500 (Tue, 05 Aug 2008) $
 * $Id: ForceCeperley.h 2945 2008-08-05 15:21:33Z jnkim $
 ***************************************************************************/

