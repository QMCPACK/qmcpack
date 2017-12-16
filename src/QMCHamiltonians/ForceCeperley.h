//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: John R. Gergely,  University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: John R. Gergely,  University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
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
  Matrix<EstimatorRealType> Sinv; // terms in fitting polynomial
  Vector<EstimatorRealType> h; // terms in fitting polynomial
  Vector<EstimatorRealType> c; // polynomial coefficients
  // container for short-range force estimator

  ParticleSet::ParticlePos_t forces_ShortRange;

  ForceCeperley(ParticleSet& ions, ParticleSet& elns);

  Return_t evaluate(ParticleSet& P);

  void InitMatrix();

  inline Return_t evaluate(ParticleSet& P, std::vector<NonLocalData>& Txy)
  {
    return evaluate(P);
  }

  void registerObservables(std::vector<observable_helper*>& h5list, hid_t gid) const
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


