//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: John R. Gergely,  University of Illinois at Urbana-Champaign
//                    Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: John R. Gergely,  University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_FORCE_BASE_HAMILTONIAN_H
#define QMCPLUSPLUS_FORCE_BASE_HAMILTONIAN_H
#include "QMCHamiltonians/ObservableHelper.h"
#include "Particle/ParticleSet.h"

namespace qmcplusplus
{
class ForceBase
{
public:
  /** cheat, need to use virtual inheriance to clean up*/
  using Real = QMCTraits::RealType;

  inline Real g(Real r)
  {
    if (r > rcut_)
      return 1.0;
    Real sum      = 0.0;
    Real r2kplusm = r;
    for (int i = 0; i < m_; i++)
      r2kplusm *= r;
    for (int k = 0; k < ck_.size(); k++)
    {
      sum += ck_[k] * r2kplusm;
      r2kplusm *= r;
    }
    return sum;
  }

  void initVarReduction(Real rcut, int m, int numFuncs);

  ForceBase(ParticleSet& ions, ParticleSet& elns);
  virtual ~ForceBase();

  void registerObservablesF(std::vector<ObservableHelper>& h5list, hdf_archive& file) const;

  void addObservablesF(QMCTraits::PropertySetType& plist);
  void addObservablesStress(QMCTraits::PropertySetType& plist);
  void setObservablesF(QMCTraits::PropertySetType& plist);
  void setObservablesStress(QMCTraits::PropertySetType& plist);
  void setParticleSetF(QMCTraits::PropertySetType& plist, int offset);
  void setParticleSetStress(QMCTraits::PropertySetType& plist, int offset);

  bool getAddIonIon() const noexcept { return add_ion_ion_; }
  void setAddIonIon(bool val) noexcept { add_ion_ion_ = val; }

  const ParticleSet::ParticlePos& getForces() const noexcept { return forces_; }
  void setForces(const ParticleSet::ParticlePos& forces);
  void setForces(Real val);

  const ParticleSet::ParticlePos& getForcesIonIon() const noexcept { return forces_ion_ion_; }
  void setForcesIonIon(const ParticleSet::ParticlePos& forces_ion_ion);

  const SymTensor<Real, OHMMS_DIM>& getStressIonIon() const noexcept { return stress_ion_ion_; }
  const SymTensor<Real, OHMMS_DIM>& getStressEE() const noexcept { return stress_ee_; }
  const SymTensor<Real, OHMMS_DIM>& getStressEI() const noexcept { return stress_ei_; }
  const SymTensor<Real, OHMMS_DIM>& getStressKin() const noexcept { return stress_kin_; }
  const SymTensor<Real, OHMMS_DIM>& getStress() const noexcept { return stress_; }

protected:
  int first_force_index_;
  int n_nuc_;
  int n_el_;
  int tries_;
  bool first_time_;
  /// Determines if ion-ion force will be added to electron-ion force in derived force estimators. If false, forces_ion_ion_=0.0.
  bool add_ion_ion_;

  ParticleSet& ions_;
  ParticleSet::ParticlePos forces_;
  ParticleSet::ParticlePos forces_ion_ion_;

  SymTensor<Real, OHMMS_DIM> stress_ion_ion_;
  SymTensor<Real, OHMMS_DIM> stress_ee_;
  SymTensor<Real, OHMMS_DIM> stress_ei_;
  SymTensor<Real, OHMMS_DIM> stress_kin_;
  SymTensor<Real, OHMMS_DIM> stress_;

  std::string prefix_;
  std::string pair_name_;

  // Data for variance reduction of Chiesa et al.
  // PRL 94, 036404 (2005)
  Real rcut_;
  int m_;
  std::vector<Real> ck_;
};

} // namespace qmcplusplus
#endif
