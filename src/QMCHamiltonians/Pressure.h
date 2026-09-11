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


#ifndef QMCPLUSPLUS_BAREPRESSURE_H
#define QMCPLUSPLUS_BAREPRESSURE_H
#include "Particle/ParticleSet.h"
#include "QMCDrivers/WalkerProperties.h"
#include "QMCHamiltonians/OperatorBase.h"


namespace qmcplusplus
{
/** @ingroup hamiltonian
 @brief Evaluate the Bare Pressure.
 P=/frac{2T+V}{d* /Omega}
 where d is the dimension of space and /Omega is the volume.
 Pressure can only be a non-Hamiltonian operator.
 In physics, it depends on the full Hamiltonian, see the above equation.
 The current implementation accesses Hamiltonian value via ParticleSet::PropertyList.
 This is a quite bad implementation choice but cannot be addressed
 until a proper ParticleSet::PropertyList gets in-place.
 **/
class Pressure : public OperatorDependsOnlyOnParticleSet
{
  using WP = WalkerProperties::Indexes;
  //     bool ZV;
  //     bool ZB;

public:
  /** constructor
   *
   * Pressure operators need to be re-evaluated during optimization.
   */
  Pressure() { update_mode_.set(OPTIMIZABLE, 1); }
  ///destructor
  ~Pressure() override {}

  std::string getClassName() const override { return "Pressure"; }

  inline Return_t evaluate(ParticleSet& P) override
  {
    const double pNorm = 1.0 / (P.getLattice().DIM * P.getLattice().Volume);
    value_             = 2.0 * P.PropertyList[WP::LOCALENERGY] - P.PropertyList[WP::LOCALPOTENTIAL];
    value_ *= 1.0 / (P.getLattice().DIM * P.getLattice().Volume);
    return 0.0;
  }

  /** implements the virtual function.
   *
   * Nothing is done but should check the mass
   */

  bool put(xmlNodePtr cur) override { return true; }

  bool get(std::ostream& os) const override
  {
    os << "Pressure";
    return true;
  }

  std::unique_ptr<OperatorBase> makeClone(ParticleSet& qp) const final { return std::make_unique<Pressure>(); }
};
} // namespace qmcplusplus
#endif
