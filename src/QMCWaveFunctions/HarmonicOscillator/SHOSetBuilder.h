//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_SHO_BASIS_BUILDER_H
#define QMCPLUSPLUS_SHO_BASIS_BUILDER_H

#include "QMCWaveFunctions/HarmonicOscillator/SHOSet.h"
#include "QMCWaveFunctions/SPOSetBuilder.h"
#include "QMCWaveFunctions/SPOSetInfo.h"

namespace qmcplusplus
{
struct SHOSetBuilder : public SPOSetBuilder
{
  //enum{DIM=OHMMS_DIM}

  ParticleSet& Ps;

  RealType length;
  RealType mass;
  RealType energy;
  PosType center;

  int nstates;
  int nmax;
  TinyVector<int, DIM> ind_dims;

  SPOSetInfoSimple<SHOState> basis_states;

  //construction/destruction
  SHOSetBuilder(ParticleSet& P, Communicate* comm);

  ~SHOSetBuilder() override;

  //reset parameters
  void reset();

  //SPOSetBuilder interface
  std::unique_ptr<SPOSet> createSPOSetFromXML(xmlNodePtr cur) override;

  std::unique_ptr<SPOSet> createSPOSet(xmlNodePtr cur, SPOSetInputInfo& input) override;

  //local functions
  void update_basis_states(int smax);
  void report(const std::string& pad = "");
};

} // namespace qmcplusplus


#endif
