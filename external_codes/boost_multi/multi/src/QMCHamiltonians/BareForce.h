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


#ifndef QMCPLUSPLUS_BARE_FORCE_HAMILTONIAN_H
#define QMCPLUSPLUS_BARE_FORCE_HAMILTONIAN_H
#include "QMCHamiltonians/ForceBase.h"
#include "QMCHamiltonians/OperatorBase.h"

namespace qmcplusplus
{
struct BareForce : public OperatorBase, public ForceBase
{
public:
  BareForce(ParticleSet& ions, ParticleSet& elns);
  std::string getClassName() const override;
  void resetTargetParticleSet(ParticleSet& P) override;

  Return_t evaluate(ParticleSet& P) override;

  void registerObservables(std::vector<ObservableHelper>& h5list, hdf_archive& file) const override;

  /** default implementation to add named values to  the property list
   * @param plist RecordNameProperty
   * @param collectables Observables that are accumulated by evaluate
   */
  void addObservables(PropertySetType& plist, BufferType& collectables) override;

  void setObservables(PropertySetType& plist) override;

  void setParticlePropertyList(PropertySetType& plist, int offset) override;

  /** Do nothing */
  bool put(xmlNodePtr cur) override;

  bool get(std::ostream& os) const override;

  std::unique_ptr<OperatorBase> makeClone(ParticleSet& qp, TrialWaveFunction& psi) final;

private:
  const int d_ei_id_;
};

} // namespace qmcplusplus
#endif
