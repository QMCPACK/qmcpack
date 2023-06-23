//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Kevin Ryczko, kryczko@uottawa.ca, University of Ottawa
//                    Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//
// File created by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_GRID_EXTERNAL_POTENTIAL_H
#define QMCPLUSPLUS_GRID_EXTERNAL_POTENTIAL_H

#include <memory>
#include "QMCHamiltonians/OperatorBase.h"
#include "einspline/bspline.h"


namespace qmcplusplus
{
/** This class allows one to read in an arbitrary external potential
  */
class GridExternalPotential : public OperatorBase
{
public:
  GridExternalPotential(ParticleSet& P);

  std::string getClassName() const override { return "GridExternalPotential"; }

  void resetTargetParticleSet(ParticleSet& P) override {}

  //standard interface functions
  bool put(xmlNodePtr cur) override;
  bool get(std::ostream& os) const override;
  std::unique_ptr<OperatorBase> makeClone(ParticleSet& P, TrialWaveFunction& psi) final;

  //functions for physical (hamiltonian component) estimator
  Return_t evaluate(ParticleSet& P) override;
  Return_t evaluate(ParticleSet& P, std::vector<NonLocalData>& Txy);

#if !defined(REMOVE_TRACEMANAGER)
  //traces interface
  void contributeParticleQuantities() override;

  void checkoutParticleQuantities(TraceManager& tm) override;

  void deleteParticleQuantities() override;

private:
  Return_t evaluate_sp(ParticleSet& P);
#endif

private:
  const ParticleSet& ps_;

  std::shared_ptr<UBspline_3d_d> spline_data_;

#if !defined(REMOVE_TRACEMANAGER)
  ///single particle trace sample array
  Array<TraceReal, 1>* v_sample_;
#endif
};
} // namespace qmcplusplus
#endif
