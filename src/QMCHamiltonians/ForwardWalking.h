//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


/** @file ForwardWalking.h
 * @brief Declarations of ForwardWalking
 */
#ifndef QMCPLUSPLUS_FORWARDWALKING_H
#define QMCPLUSPLUS_FORWARDWALKING_H
#include "QMCHamiltonians/OperatorBase.h"

namespace qmcplusplus
{
class QMCHamiltonian;

class ForwardWalking : public OperatorBase
{
public:
  /** constructor
   */
  ForwardWalking();

  ///destructor
  ~ForwardWalking() override;

  std::string getClassName() const override { return "ForwardWalking"; }

  void resetTargetParticleSet(ParticleSet& P) override {}

  Return_t rejectedMove(ParticleSet& P) override;

  Return_t calculate(ParticleSet& P);

  Return_t evaluate(ParticleSet& P) override;

  bool put(xmlNodePtr cur) override;

  ///rename it to avoid conflicts with put
  bool putSpecial(xmlNodePtr cur, QMCHamiltonian& h, ParticleSet& P);

  bool get(std::ostream& os) const override;

  std::unique_ptr<OperatorBase> makeClone(ParticleSet& qp, TrialWaveFunction& psi) final;

  void addObservables(PropertySetType& plist);

  void addObservables(PropertySetType& plist, BufferType& collectables) override;

  void setObservables(PropertySetType& plist) override;

  void setParticlePropertyList(PropertySetType& plist, int offset) override;

private:
  std::vector<int> h_ids_;
  std::vector<int> p_ids_;
  std::vector<std::vector<int>> walker_lengths_;
  std::vector<RealType> values_;
  std::vector<std::string> names_;
  int nobservables_;
  int nvalues_;
  int first_hamiltonian_;
};
} // namespace qmcplusplus
#endif
