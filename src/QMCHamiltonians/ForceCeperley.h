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
#include "QMCHamiltonians/OperatorBase.h"
#include "LongRange/LRCoulombSingleton.h"
#include "Numerics/OneDimGridBase.h"
#include "Numerics/OneDimGridFunctor.h"
#include "Numerics/OneDimCubicSpline.h"

namespace qmcplusplus
{
struct ForceCeperley : public OperatorBase, public ForceBase
{
private:
  const int d_aa_ID;
  const int d_ei_ID;

public:
  double Rcut;                   // parameter: radial distance within which estimator is used
  int m_exp;                     // parameter: exponent in polynomial fit
  int N_basis;                   // parameter: size of polynomial basis set
  Matrix<FullPrecRealType> Sinv; // terms in fitting polynomial
  Vector<FullPrecRealType> h;    // terms in fitting polynomial
  Vector<FullPrecRealType> c;    // polynomial coefficients
  // container for short-range force estimator

  ForceCeperley(ParticleSet& ions, ParticleSet& elns);

  std::string getClassName() const override { return "ForceCeperley"; }

  Return_t evaluate(ParticleSet& P) override;

  void InitMatrix();

  void registerObservables(std::vector<ObservableHelper>& h5list, hdf_archive& file) const override
  {
    registerObservablesF(h5list, file);
  }

  void addObservables(PropertySetType& plist, BufferType& collectables) override { addObservablesF(plist); }

  void setObservables(PropertySetType& plist) override { setObservablesF(plist); }

  void resetTargetParticleSet(ParticleSet& P) override {}

  // Compute ion-ion forces at construction to include in the total forces
  void evaluate_IonIon(ParticleSet::ParticlePos& forces) const;

  void setParticlePropertyList(PropertySetType& plist, int offset) override { setParticleSetF(plist, offset); }
  std::unique_ptr<OperatorBase> makeClone(ParticleSet& qp, TrialWaveFunction& psi) final;

  bool put(xmlNodePtr cur) override;

  bool get(std::ostream& os) const override
  {
    os << "Ceperley Force Estimator Hamiltonian: " << pair_name_;
    return true;
  }
};

} // namespace qmcplusplus
#endif
