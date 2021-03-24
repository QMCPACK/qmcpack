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

#include "Configuration.h"
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
#include "Utilities/StlPrettyPrint.hpp"

namespace qmcplusplus
{
// The missing unit tests for TrialWaveFunction.
// Right now only uses the "primary" DiracDeterminant type.
#if defined(ENABLE_CUDA) && defined(ENABLE_OFFLOAD)
using VT       = QMCTraits::ValueType;
using FPVT     = QMCTraits::QTFull::ValueType;
using DiracDet = DiracDeterminantBatched<MatrixDelayedUpdateCUDA<VT, FPVT>>;
#else
#ifdef ENABLE_CUDA
using DiracDet = DiracDeterminant<DelayedUpdateCUDA<QMCTraits::ValueType, QMCTraits::QTFull::ValueType>>;
#else
using DiracDet = DiracDeterminant<DelayedUpdate<QMCTraits::ValueType, QMCTraits::QTFull::ValueType>>;
#endif
#endif

using LogValueType = TrialWaveFunction::LogValueType;
using PsiValueType = TrialWaveFunction::PsiValueType;
using GradType     = TrialWaveFunction::GradType;
using PosType      = QMCTraits::PosType;
using RealType     = QMCTraits::RealType;
using ValueType    = QMCTraits::ValueType;

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

class MinimalTrialWaveFunction
{
public:
  MinimalTrialWaveFunction()
  {
    Communicate* c = OHMMS::Controller;

    ions = std::make_unique<ParticleSet>();
    elec = std::make_unique<ParticleSet>();
    ParticleSet& ions_(*ions);
    ParticleSet& elec_(*elec);

    ions_.setName("ion");
    ions_.create(2);
    ions_.R[0] = {0.0, 0.0, 0.0};
    ions_.R[1] = {1.68658058, 1.68658058, 1.68658058};
    ions_.update();

    elec_.setName("elec");
    std::vector<int> ud(2);
    ud[0] = ud[1] = 2;
    elec_.create(ud);
    elec_.R[0] = {0.0, 0.0, 0.0};
    elec_.R[1] = {0.0, 1.0, 1.0};
    elec_.R[2] = {1.0, 1.0, 0.0};
    elec_.R[3] = {1.0, 0.0, 1.0};

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

    // make a ParticleSet Clone
    ParticleSet elec_clone(elec_);

    //diamondC_1x1x1
    const char* spo_xml = "<tmp> \
<determinantset type=\"einspline\" href=\"diamondC_1x1x1.pwscf.h5\" tilematrix=\"1 0 0 0 1 0 0 0 1\" twistnum=\"0\" source=\"ion\" meshfactor=\"1.0\" precision=\"float\" size=\"2\"/> \
</tmp> \
";

    Libxml2Document doc;
    bool okay = doc.parseFromString(spo_xml);
    REQUIRE(okay);

    xmlNodePtr spo_root = doc.getRoot();
    xmlNodePtr ein1     = xmlFirstElementChild(spo_root);

    // EinsplineSetBuilder requires the pool and holds on to it.
    UPtr<ParticleSet> elec_uptr = std::make_unique<ParticleSet>(elec_);
    UPtr<ParticleSet> ions_uptr = std::make_unique<ParticleSet>(ions_);
    ParticleSetPool ptcl{c};
    ptcl.addParticleSet(std::move(elec_uptr));
    ptcl.addParticleSet(std::move(ions_uptr));

    EinsplineSetBuilder einSet(elec_, ptcl.getPool(), c, ein1);
    spo = einSet.createSPOSetFromXML(ein1);
    REQUIRE(spo != nullptr);

    // Make DiracDets
    det_up = new DiracDet(std::unique_ptr<SPOSet>(spo->makeClone()));
    det_up->set(0, 2);
    det_dn = new DiracDet(std::unique_ptr<SPOSet>(spo->makeClone()));
    det_dn->set(2, 2);

    // Make SlaterDet
    slater_det = new SlaterDet(elec_);
    slater_det->add(det_up, 0);
    slater_det->add(det_dn, 1);

    // Make TWFC with just slater det.
    psi = std::make_unique<TrialWaveFunction>();
    psi->addComponent(slater_det);
  }

  SlaterDet* slater_det;
  DiracDet* det_up;
  DiracDet* det_dn;
  UPtr<ParticleSet> ions;
  UPtr<ParticleSet> elec;
  UPtr<TrialWaveFunction> psi;

  SPOSet* spo;
};

// Doesn't support legacy CUDA
#if !defined(QMC_CUDA)

TEST_CASE("TrialWaveFunction_evaluateLog", "[wavefunction][fermion]")
{
  MinimalTrialWaveFunction mtw;

  ResourceCollection res_col("test");
  TrialWaveFunction& psi = *(mtw.psi);
  psi.createResource(res_col);
  psi.acquireResource(res_col);

  auto& elec = *(mtw.elec);
  elec.update();

  ParticleSet elec_clone(elec);

  double logpsi = psi.evaluateLog(elec);
  CHECK(logpsi == psi.getLogPsi());
  CHECK(std::complex<RealType>(psi.getLogPsi(), psi.getPhase()) ==
        LogComplexApprox(std::complex<RealType>(-3.75034, 9.42478)));
}

TEST_CASE("TrialWaveFunction_mw_evaluateLog", "[wavefunction][fermion]")
{
  MinimalTrialWaveFunction mtw;

  ResourceCollection res_col("test");
  TrialWaveFunction& psi = *(mtw.psi);
  psi.createResource(res_col);
  psi.acquireResource(res_col);

  auto& elec = *(mtw.elec);
  elec.update();

  ParticleSet elec_clone(elec);

  UPtr<TrialWaveFunction> psi_clone_ptr(psi.makeClone(elec_clone));
  TrialWaveFunction& psi_clone = *psi_clone_ptr;

  RefVectorWithLeader<ParticleSet> p_ref_list(elec, {elec, elec_clone});
  RefVectorWithLeader<TrialWaveFunction> wf_ref_list(psi, {psi, psi_clone});

  ParticleSet::mw_update(p_ref_list);
  TrialWaveFunction::mw_evaluateLog(wf_ref_list, p_ref_list);
  CHECK(std::complex<RealType>(wf_ref_list[0].getLogPsi(), wf_ref_list[0].getPhase()) ==
        LogComplexApprox(std::complex<RealType>(-3.75034, 9.42478)));
  CHECK(std::complex<RealType>(wf_ref_list[1].getLogPsi(), wf_ref_list[1].getPhase()) ==
        LogComplexApprox(std::complex<RealType>(-3.75034, 9.42478)));
  res_col.rewind();
  psi.releaseResource(res_col);
  res_col.rewind();
}

TEST_CASE("TrialWaveFunction_evaluateGrad", "[wavefunction][fermion]")
{
  MinimalTrialWaveFunction mtw;

  ResourceCollection res_col("test");
  TrialWaveFunction& psi = *(mtw.psi);
  psi.createResource(res_col);
  psi.acquireResource(res_col);

  auto& elec = *(mtw.elec);
  elec.update();

  double logpsi = psi.evaluateLog(elec);

  using PosType   = QMCTraits::PosType;
  using RealType  = QMCTraits::RealType;
  using ValueType = QMCTraits::ValueType;
  PosType delta(0.1, 0.1, 0.2);
  int moved_elec_id = 0;
  elec.makeMove(moved_elec_id, delta);

  GradType grad = psi.evalGrad(elec, 0);
  std::cout << grad << '\n';
}

TEST_CASE("TrialWaveFunction_mw_evaluateGrad", "[wavefunction][fermion]")
{
  MinimalTrialWaveFunction mtw;
  ResourceCollection res_col("test");
  TrialWaveFunction& psi = *(mtw.psi);
  psi.createResource(res_col);
  psi.acquireResource(res_col);

  auto& elec = *(mtw.elec);
  elec.update();

  ParticleSet elec_clone(elec);

  UPtr<TrialWaveFunction> psi_clone_ptr(psi.makeClone(elec_clone));
  TrialWaveFunction& psi_clone = *psi_clone_ptr;

  RefVectorWithLeader<ParticleSet> p_ref_list(elec, {elec, elec_clone});
  RefVectorWithLeader<TrialWaveFunction> wf_ref_list(psi, {psi, psi_clone});
  int nw = p_ref_list.size();
  ParticleSet::mw_update(p_ref_list);
  TrialWaveFunction::mw_evaluateLog(wf_ref_list, p_ref_list);

  using PosType   = QMCTraits::PosType;
  using RealType  = QMCTraits::RealType;
  using ValueType = QMCTraits::ValueType;
  std::vector<PosType> deltas{{0.1, 0.1, 0.2}, {0.05, 0.15, 0.15}};
  int moved_elec_id = 0;
  ParticleSet::mw_makeMove(p_ref_list, moved_elec_id, deltas);
  std::vector<GradType> grads_now(nw);
  TrialWaveFunction::mw_evalGrad(wf_ref_list, p_ref_list, moved_elec_id, grads_now);
  std::cout << grads_now << '\n';
}

#endif

} // namespace qmcplusplus
