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
#include "QMCHamiltonians/OperatorBase.h"

namespace qmcplusplus
{
struct ForceBase
{
  /** cheat, need to use virtual inheriance to clean up*/
  using real_type = QMCTraits::RealType;

  int FirstForceIndex;
  int Nnuc;
  int Nel;
  int tries;
  bool FirstTime;
  bool addionion;

  ParticleSet& Ions;
  ParticleSet::ParticlePos forces;
  ParticleSet::ParticlePos forces_IonIon;
  SymTensor<real_type, OHMMS_DIM> stress_IonIon;
  SymTensor<real_type, OHMMS_DIM> stress_ee, stress_ei, stress_kin;
  SymTensor<real_type, OHMMS_DIM> stress;

  std::string prefix;
  std::string pairName;

  // Data for variance reduction of Chiesa et al.
  // PRL 94, 036404 (2005)
  real_type Rcut;
  int m;
  std::vector<real_type> ck;
  inline real_type g(real_type r)
  {
    if (r > Rcut)
      return 1.0;
    real_type sum      = 0.0;
    real_type r2kplusm = r;
    for (int i = 0; i < m; i++)
      r2kplusm *= r;
    for (int k = 0; k < ck.size(); k++)
    {
      sum += ck[k] * r2kplusm;
      r2kplusm *= r;
    }
    return sum;
  }


  void InitVarReduction(real_type Rcut, int m, int numFuncs);


  ForceBase(ParticleSet& ions, ParticleSet& elns);
  virtual ~ForceBase() {}

  void registerObservablesF(std::vector<ObservableHelper>& h5list, hid_t gid) const;

  void addObservablesF(QMCTraits::PropertySetType& plist);
  void addObservablesStress(QMCTraits::PropertySetType& plist);
  void setObservablesF(QMCTraits::PropertySetType& plist);
  void setObservablesStress(QMCTraits::PropertySetType& plist);
  void setParticleSetF(QMCTraits::PropertySetType& plist, int offset);
  void setParticleSetStress(QMCTraits::PropertySetType& plist, int offset);
};

struct BareForce : public OperatorBase, public ForceBase
{
private:
  const int d_ei_ID;

public:
  BareForce(ParticleSet& ions, ParticleSet& elns);
  void resetTargetParticleSet(ParticleSet& P) override;

  Return_t evaluate(ParticleSet& P) override;

  void registerObservables(std::vector<ObservableHelper>& h5list, hid_t gid) const override
  {
    registerObservablesF(h5list, gid);
  }

  /** default implementation to add named values to  the property list
   * @param plist RecordNameProperty
   * @param collectables Observables that are accumulated by evaluate
   */
  void addObservables(PropertySetType& plist, BufferType& collectables) override;

  void setObservables(PropertySetType& plist) override { setObservablesF(plist); }

  void setParticlePropertyList(PropertySetType& plist, int offset) override { setParticleSetF(plist, offset); }

  /** Do nothing */
  bool put(xmlNodePtr cur) override;

  bool get(std::ostream& os) const override
  {
    os << "Force Base Hamiltonian: " << pairName;
    return true;
  }

  std::unique_ptr<OperatorBase> makeClone(ParticleSet& qp, TrialWaveFunction& psi) final;
};

} // namespace qmcplusplus
#endif
