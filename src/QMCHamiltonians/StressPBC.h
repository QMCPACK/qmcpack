//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers
//
// File developed by: Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Yubo Paul Yang, yyang173@illinois.edu, University of Illinois Urbana-Champaign
//
// File created by: Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_STRESSPBC_HAMILTONIAN_H
#define QMCPLUSPLUS_STRESSPBC_HAMILTONIAN_H
#include "QMCHamiltonians/ForceBase.h"
#include "QMCHamiltonians/OperatorBase.h"
#include "LongRange/LRCoulombSingleton.h"
#include "OhmmsPETE/SymTensor.h"

namespace qmcplusplus
{
struct StressPBC : public OperatorBase, public ForceBase
{
  using LRHandlerType = LRCoulombSingleton::LRHandlerType;

  SymTensor<RealType, OHMMS_DIM> stress_ee_const;
  SymTensor<RealType, OHMMS_DIM> stress_eI_const;

  // need wave function for kinetic stress
  TrialWaveFunction& Psi;

  ///source particle set
  ParticleSet& PtclTarg;
  ParticleSet& PtclA;
  ///long-range Handler
  std::unique_ptr<LRHandlerType> AA;
  ///locator of the distance table
  const int ei_table_index;
  /// e-e table ID
  const int ee_table_index;
  const int ii_table_index;
  ///number of species of A particle set
  int NumSpeciesA;
  ///number of species of B particle set
  int NumSpeciesB;
  ///number of particles of A (classical, e.g. ions)
  int NptclA;
  ///number of particles of B (quantum, e.g. electrons)
  int NptclB;

  ///number of particles per species of A
  std::vector<int> NofSpeciesA;
  ///number of particles per species of B
  std::vector<int> NofSpeciesB;
  ///Zat[iat] charge for the iat-th particle of A
  std::vector<RealType> Zat;
  ///Qat[iat] charge for the iat-th particle of B
  std::vector<RealType> Qat;
  ///Zspec[spec] charge for the spec-th species of A
  std::vector<RealType> Zspec;
  ///Qspec[spec] charge for the spec-th species of B
  std::vector<RealType> Qspec;
  //Constructor
  bool firstTimeStress;
  StressPBC(ParticleSet& ions, ParticleSet& elns, TrialWaveFunction& Psi);

  std::string getClassName() const override { return "StressPBC"; }

  Return_t evaluate(ParticleSet& P) override;

  void initBreakup(ParticleSet& P);

  SymTensor<RealType, OHMMS_DIM> evaluateLR_AB(ParticleSet& P);
  SymTensor<RealType, OHMMS_DIM> evaluateSR_AB(ParticleSet& P_target);
  SymTensor<RealType, OHMMS_DIM> evaluateSR_AA(ParticleSet& P, int itabSelf);
  SymTensor<RealType, OHMMS_DIM> evaluateLR_AA(ParticleSet& P);
  SymTensor<RealType, OHMMS_DIM> evalConsts_AB();
  SymTensor<RealType, OHMMS_DIM> evalConsts_AA(ParticleSet& P);

  SymTensor<RealType, OHMMS_DIM> evaluateKineticSymTensor(ParticleSet& P);

  void registerObservables(std::vector<ObservableHelper>& h5list, hdf_archive& file) const override
  {
    registerObservablesF(h5list, file);
  }

  void addObservables(PropertySetType& plist, BufferType& collectables) override { addObservablesStress(plist); }

  void setObservables(PropertySetType& plist) override { setObservablesStress(plist); }

  void resetTargetParticleSet(ParticleSet& P) override {}

  void setParticlePropertyList(PropertySetType& plist, int offset) override { setParticleSetStress(plist, offset); }
  std::unique_ptr<OperatorBase> makeClone(ParticleSet& qp, TrialWaveFunction& psi) final;
  bool put(xmlNodePtr cur) override;

  bool get(std::ostream& os) const override
  {
    os << "Ceperley Force Estimator Hamiltonian: " << pair_name_;
    return true;
  }

  void CalculateIonIonStress()
  {
    stress_ion_ion_ = evaluateSR_AA(PtclA, ii_table_index) + evaluateLR_AA(PtclA) + evalConsts_AA(PtclA);
    stress_eI_const += evalConsts_AB();
    stress_ee_const += evalConsts_AA(PtclTarg);
  }
};

} // namespace qmcplusplus
#endif
