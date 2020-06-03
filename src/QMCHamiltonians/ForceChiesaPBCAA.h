//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_FORCE_CHIESA_HAMILTONIAN_H
#define QMCPLUSPLUS_FORCE_CHIESA_HAMILTONIAN_H
#include "QMCHamiltonians/ForceBase.h"
#include "LongRange/LRCoulombSingleton.h"
#include "Numerics/OneDimGridBase.h"
#include "Numerics/OneDimGridFunctor.h"
#include "Numerics/OneDimCubicSpline.h"

namespace qmcplusplus
{
struct ForceChiesaPBCAA : public OperatorBase, public ForceBase
{
  typedef LRCoulombSingleton::LRHandlerType LRHandlerType;
  typedef LRCoulombSingleton::GridType GridType;
  typedef LRCoulombSingleton::RadFunctorType RadFunctorType;

  RealType Rcut;         // parameter: radial distance within which estimator is used
  int m_exp;             // parameter: exponent in polynomial fit
  int N_basis;           // parameter: size of polynomial basis set
  Matrix<RealType> Sinv; // terms in fitting polynomial
  Vector<RealType> h;    // terms in fitting polynomial
  Vector<RealType> c;    // polynomial coefficients
  // container for short-range force estimator

  ///source particle set
  ParticleSet& PtclA;
  ///long-range Handler
  LRHandlerType* dAB;
  ///number of species of A particle set
  int NumSpeciesA;
  ///number of species of B particle set
  int NumSpeciesB;
  ///number of particles of A
  int NptclA;
  ///number of particles of B
  int NptclB;

  ///Zat[iat] charge for the iat-th particle of A
  std::vector<RealType> Zat;
  ///Qat[iat] charge for the iat-th particle of B
  std::vector<RealType> Qat;
  ///Zspec[spec] charge for the spec-th species of A
  std::vector<RealType> Zspec;
  ///Qspec[spec] charge for the spec-th species of B
  std::vector<RealType> Qspec;

  bool first_time;

  ForceChiesaPBCAA(ParticleSet& ions, ParticleSet& elns, bool firsttime = true);

  Return_t evaluate(ParticleSet& P);

  void InitMatrix();
  void initBreakup(ParticleSet& P);

  void evaluateLR(ParticleSet&);
  void evaluateSR(ParticleSet&);
  void evaluateSR_AA();
  void evaluateLR_AA();

  Return_t g_filter(RealType r);

  void registerObservables(std::vector<observable_helper*>& h5list, hid_t gid) const
  {
    registerObservablesF(h5list, gid);
  }

  void addObservables(PropertySetType& plist, BufferType& collectables);


  void setObservables(PropertySetType& plist)
  {
    OperatorBase::setObservables(plist);
    setObservablesF(plist);
  }

  void setParticlePropertyList(PropertySetType& plist, int offset)
  {
    OperatorBase::setParticlePropertyList(plist, offset);
    setParticleSetF(plist, offset);
  }


  void resetTargetParticleSet(ParticleSet& P);


  OperatorBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi);

  bool put(xmlNodePtr cur);

  bool get(std::ostream& os) const
  {
    os << "Ceperley Force Estimator Hamiltonian: " << pairName;
    return true;
  }

  // for testing only
  int getDistanceTableAAID() const { return d_aa_ID; }

private:
  // AA table ID
  const int d_aa_ID;
  const int d_ei_ID;
};

} // namespace qmcplusplus
#endif
