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
#include "QMCHamiltonians/OperatorBase.h"
#include "LongRange/LRCoulombSingleton.h"
#include "Numerics/OneDimGridBase.h"
#include "Numerics/OneDimGridFunctor.h"
#include "Numerics/OneDimCubicSpline.h"

namespace qmcplusplus
{
struct ForceChiesaPBCAA : public OperatorBase, public ForceBase
{
  using LRHandlerType  = LRCoulombSingleton::LRHandlerType;
  using GridType       = LRCoulombSingleton::GridType;
  using RadFunctorType = LRCoulombSingleton::RadFunctorType;

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
  std::unique_ptr<LRHandlerType> dAB;
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

  std::string getClassName() const override { return "ForceChiesaPBCAA"; }

  Return_t evaluate(ParticleSet& P) override;

  void InitMatrix();
  void initBreakup(ParticleSet& P);

  void evaluateLR(ParticleSet&);
  void evaluateSR(ParticleSet&);
  void evaluateSR_AA();
  void evaluateLR_AA();

  Return_t g_filter(RealType r);

  void registerObservables(std::vector<ObservableHelper>& h5list, hdf_archive& file) const override
  {
    registerObservablesF(h5list, file);
  }

  void addObservables(PropertySetType& plist, BufferType& collectables) override;


  void setObservables(PropertySetType& plist) override
  {
    OperatorBase::setObservables(plist);
    setObservablesF(plist);
  }

  void setParticlePropertyList(PropertySetType& plist, int offset) override
  {
    OperatorBase::setParticlePropertyList(plist, offset);
    setParticleSetF(plist, offset);
  }


  void resetTargetParticleSet(ParticleSet& P) override;


  std::unique_ptr<OperatorBase> makeClone(ParticleSet& qp, TrialWaveFunction& psi) final;

  bool put(xmlNodePtr cur) override;

  bool get(std::ostream& os) const override
  {
    os << "Ceperley Force Estimator Hamiltonian: " << pair_name_;
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
