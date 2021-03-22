//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#include <catch.hpp>

#include "OhmmsData/Libxml2Doc.h"
#include "Particle/ParticleSet.h"
#include "Particle/ParticleSetPool.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "QMCWaveFunctions/EinsplineSetBuilder.h"
#include "QMCWaveFunctions/Fermion/DiracDeterminant.h"
#include "QMCWaveFunctions/Fermion/DiracDeterminantBatched.h"
#include "QMCWaveFunctions/Fermion/SlaterDet.h"
#include "QMCWaveFunctions/Jastrow/RadialJastrowBuilder.h"
#include "ResourceCollection.h"
#include "FakeSPO.h"

namespace qmcplusplus
{
// The missing unit tests for TrialWaveFunction.

#if defined(ENABLE_CUDA) && defined(ENABLE_OFFLOAD)
using VT   = QMCTraits::ValueType;
using FPVT = QMCTraits::QTFull::ValueType;
using DiracDet = DiracDeterminantBatched<MatrixDelayedUpdateCUDA<VT, FPVT>>;
#else
#ifdef ENABLE_CUDA
using DiracDet = DiracDeterminant<DelayedUpdateCUDA<QMCTraits::ValueType, QMCTraits::QTFull::ValueType>>;
#else
using DiracDet = DiracDeterminant<DelayedUpdate<QMCTraits::ValueType, QMCTraits::QTFull::ValueType>>;
#endif
#endif

TEST_CASE("TrialWaveFunction_resource_loop", "[wavefunction][fermion]")
{
  // It should take much less than this to setup a "valid Trialwavefunction" to get a
  // reasonable unit test and illustration of usage of
  // ::[create|acquire|release|]Resource in trialwavefunction.
  Communicate* c = OHMMS::Controller;

  auto ions_uptr = std::make_unique<ParticleSet>();
  auto elec_uptr = std::make_unique<ParticleSet>();
  ParticleSet& ions_(*ions_uptr);
  ParticleSet& elec_(*elec_uptr);

  ions_.setName("ion");
  ions_.create(2);
  ions_.R[0] = {0.0, 0.0, 0.0};
  ions_.R[1] = {1.68658058, 1.68658058, 1.68658058};
  ions_.update();


  elec_.setName("elec");
  std::vector<int> ud(2);
  ud[0] = ud[1] = 3;
  elec_.create(ud);
  elec_.R[0] = {0.0, 0.0, 0.0};
  elec_.R[1] = {0.0, 1.0, 1.0};
  elec_.R[2] = {1.0, 1.0, 0.0};
  elec_.R[3] = {1.0, 0.0, 1.0};
  elec_.R[4] = {0.5, 1.0, 0.0};
  elec_.R[5] = {1.0, 0.5, 1.0};

  // diamondC_1x1x1
  elec_.Lattice.R(0, 0)   = 3.37316115;
  elec_.Lattice.R(0, 1)   = 3.37316115;
  elec_.Lattice.R(0, 2)   = 0.0;
  elec_.Lattice.R(1, 0)   = 0.0;
  elec_.Lattice.R(1, 1)   = 3.37316115;
  elec_.Lattice.R(1, 2)   = 3.37316115;
  elec_.Lattice.R(2, 0)   = 3.37316115;
  elec_.Lattice.R(2, 1)   = 0.0;
  elec_.Lattice.R(2, 2)   = 3.37316115;
  elec_.Lattice.BoxBConds = {1, 1, 1};
  elec_.Lattice.reset();

  SpeciesSet& tspecies         = elec_.getSpeciesSet();
  int upIdx                    = tspecies.addSpecies("u");
  int downIdx                  = tspecies.addSpecies("d");
  int chargeIdx                = tspecies.addAttribute("charge");
  tspecies(chargeIdx, upIdx)   = -1;
  tspecies(chargeIdx, downIdx) = -1;

  elec_.addTable(ions_);
  elec_.resetGroups();

  ParticleSetPool ptcl{c};
  ptcl.addParticleSet(std::move(elec_uptr));
  ptcl.addParticleSet(std::move(ions_uptr));

  // make a ParticleSet Clone
  ParticleSet elec_clone(elec_);

  auto spo_init = std::make_unique<FakeSPO>();
  spo_init->setOrbitalSetSize(3);
  DiracDet ddb(std::move(spo_init));
  auto spo = dynamic_cast<FakeSPO*>(ddb.getPhi());

  auto* det_up = new DiracDet(std::unique_ptr<SPOSet>(spo->makeClone()));
  det_up->set(0, 3);
  auto* det_dn = new DiracDet(std::unique_ptr<SPOSet>(spo->makeClone()));
  det_dn->set(3, 3);

  auto* slater_det = new SlaterDet(elec_);
  slater_det->add(det_up, 0);
  slater_det->add(det_dn, 1);

  TrialWaveFunction psi;
  psi.addComponent(slater_det);

  ResourceCollection res_col("test");
  psi.createResource(res_col);
  psi.acquireResource(res_col);
  res_col.rewind();
  psi.releaseResource(res_col);

  // hand the resources back and forth;

  std::unique_ptr<TrialWaveFunction> psi_clone(psi.makeClone(elec_clone));
  res_col.rewind();
  psi_clone->acquireResource(res_col);
  res_col.rewind();  
  psi_clone->releaseResource(res_col);
  res_col.rewind();
  psi.acquireResource(res_col);
  
}

}
