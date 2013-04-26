//////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Jeongnim Kim and Bryan Clark
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
#ifndef QMCPLUSPLUS_ONEBODYDENSITYMATRIX_HAMILTONIAN_H
#define QMCPLUSPLUS_ONEBODYDENSITYMATRIX_HAMILTONIAN_H
#include "QMCHamiltonians/QMCHamiltonianBase.h"

namespace qmcplusplus
{

class OneBodyDensityMatrix: public QMCHamiltonianBase
{
public:
  OneBodyDensityMatrix(ParticleSet& ions, ParticleSet& elns);

  OneBodyDensityMatrix(ParticleSet& elns);

  void resetTargetParticleSet(ParticleSet& P);

  Return_t evaluate(ParticleSet& P);

  inline Return_t evaluate(ParticleSet& P, vector<NonLocalData>& Txy)
  {
    return evaluate(P);
  }

  void addObservables(PropertySetType& plist);
  void setObservables(PropertySetType& plist);
  void setParticlePropertyList(PropertySetType& plist, int offset);
  bool put(xmlNodePtr cur);
  bool get(std::ostream& os) const;
  QMCHamiltonianBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi);

private:
  /// maximum distance
  RealType Dmax;
  /// bin size
  RealType Delta;
  /// one of bin size
  RealType DeltaInv;
  ///trial wavefunction
  TrialWaveFunction& Psi;
  ///data
  Vector<RealType> gofr;
  /** resize the internal data
   *
   * The argument list is not completed
   */
  void resize();
};

}
#endif

/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 2945 $   $Date: 2008-08-05 10:21:33 -0500 (Tue, 05 Aug 2008) $
 * $Id: ForceBase.h 2945 2008-08-05 15:21:33Z jnkim $
 ***************************************************************************/
