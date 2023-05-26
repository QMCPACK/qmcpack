//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Peter W. Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "CoulombPBCAA.h"
#include <numeric>
#include "EwaldRef.h"
#include "Particle/DistanceTable.h"
#include "Utilities/ProgressReportEngine.h"
#include <ResourceCollection.h>
#include <Message/UniformCommunicateError.h>
#include "Numerics/OneDimCubicSplineLinearGrid.h"

namespace qmcplusplus
{
struct CoulombPBCAA::CoulombPBCAAMultiWalkerResource : public Resource
{
  CoulombPBCAAMultiWalkerResource() : Resource("CoulombPBCAA") {}

  std::unique_ptr<Resource> makeClone() const override
  {
    return std::make_unique<CoulombPBCAAMultiWalkerResource>(*this);
  }

  Vector<CoulombPBCAA::Return_t, OffloadPinnedAllocator<CoulombPBCAA::Return_t>> values_offload;

  /// a walkers worth of per particle coulomb AA potential values
  Vector<RealType> v_sample;

  /// constant values per particle for coulomb AA potential
  Vector<RealType> pp_consts;
};

CoulombPBCAA::CoulombPBCAA(ParticleSet& ref, bool active, bool computeForces, bool use_offload)
    : ForceBase(ref, ref),
      is_active(active),
      FirstTime(true),
      myConst(0.0),
      ComputeForces(computeForces),
      quasi2d(LRCoulombSingleton::this_lr_type == LRCoulombSingleton::QUASI2D),
      Ps(ref),
      use_offload_(active && !computeForces && use_offload),
      d_aa_ID(ref.addTable(ref, use_offload_ ? DTModes::ALL_OFF : DTModes::NEED_FULL_TABLE_ON_HOST_AFTER_DONEPBYP)),
      evalLR_timer_(createGlobalTimer("CoulombPBCAA::LongRange", timer_level_fine)),
      evalSR_timer_(createGlobalTimer("CoulombPBCAA::ShortRange", timer_level_fine)),
      offload_timer_(createGlobalTimer("CoulombPBCAA::offload", timer_level_fine))
{
  if (use_offload_)
    assert(ref.getCoordinates().getKind() == DynamicCoordinateKind::DC_POS_OFFLOAD);

  ReportEngine PRE("CoulombPBCAA", "CoulombPBCAA");
  setEnergyDomain(POTENTIAL);
  twoBodyQuantumDomain(ref);
  PtclRefName = ref.getDistTable(d_aa_ID).getName();
  if (ComputeForces || quasi2d)
  {
    ref.turnOnPerParticleSK();
  }
  initBreakup(ref);
  if (ComputeForces)
  {
    updateSource(ref);
  }

  if (!is_active)
  {
    ref.update();
    updateSource(ref);

    ewaldref::RealMat A;
    ewaldref::PosArray R;
    ewaldref::ChargeArray Q;

    A = Ps.getLattice().R;

    R.resize(NumCenters);
    Q.resize(NumCenters);
    for (int i = 0; i < NumCenters; ++i)
    {
      R[i] = Ps.R[i];
      Q[i] = Zat[i];
    }

    RealType Vii_ref        = ewaldref::ewaldEnergy(A, R, Q);
    RealType Vdiff_per_atom = std::abs(value_ - Vii_ref) / NumCenters;
    app_log() << "Checking ion-ion Ewald energy against reference..." << std::endl;
    if (Vdiff_per_atom > Ps.getLattice().LR_tol)
    {
      std::ostringstream msg;
      msg << std::setprecision(14);
      msg << "in ion-ion Ewald energy exceeds " << Ps.getLattice().LR_tol << " Ha/atom tolerance." << std::endl;
      msg << std::endl;
      msg << "  Reference ion-ion energy: " << Vii_ref << std::endl;
      msg << "  QMCPACK   ion-ion energy: " << value_ << std::endl;
      msg << "            ion-ion diff  : " << value_ - Vii_ref << std::endl;
      msg << "            diff/atom     : " << (value_ - Vii_ref) / NumCenters << std::endl;
      msg << "            tolerance     : " << Ps.getLattice().LR_tol << std::endl;
      msg << std::endl;
      msg << "Please try increasing the LR_dim_cutoff parameter in the <simulationcell/>" << std::endl;
      msg << "input.  Alternatively, the tolerance can be increased by setting the" << std::endl;
      msg << "LR_tol parameter in <simulationcell/> to a value greater than " << Ps.getLattice().LR_tol << ". "
          << std::endl;
      msg << "If you increase the tolerance, please perform careful checks of energy" << std::endl;
      msg << "differences to ensure this error is controlled for your application." << std::endl;
      msg << std::endl;

      throw UniformCommunicateError(msg.str());
    }
    else
    {
      app_log() << "  Check passed." << std::endl;
    }
  }
  prefix_ = "F_AA";
  app_log() << "  Maximum K shell " << AA->MaxKshell << std::endl;
  app_log() << "  Number of k vectors " << AA->Fk.size() << std::endl;
  app_log() << "  Fixed Coulomb potential for " << ref.getName();
  app_log() << "\n    e-e Madelung Const. =" << std::setprecision(8) << madelung_constant_
            << "\n    Vtot     =" << value_ << std::endl;
}

CoulombPBCAA::~CoulombPBCAA() = default;

void CoulombPBCAA::addObservables(PropertySetType& plist, BufferType& collectables)
{
  addValue(plist);
  if (ComputeForces)
    addObservablesF(plist);
}

void CoulombPBCAA::updateSource(ParticleSet& s)
{
  mRealType eL(0.0), eS(0.0);
  if (ComputeForces)
  {
    forces_ = 0.0;
    eS      = evalSRwithForces(s);
    eL      = evalLRwithForces(s);
  }
  else
  {
    eL = evalLR(s);
    eS = evalSR(s);
  }
  new_value_ = value_ = eL + eS + myConst;
}

void CoulombPBCAA::resetTargetParticleSet(ParticleSet& P)
{
  if (is_active)
  {
    PtclRefName = P.getDistTable(d_aa_ID).getName();
    AA->resetTargetParticleSet(P);
  }
}


#if !defined(REMOVE_TRACEMANAGER)
void CoulombPBCAA::contributeParticleQuantities() { request_.contribute_array(name_); }

void CoulombPBCAA::checkoutParticleQuantities(TraceManager& tm)
{
  streaming_particles_ = request_.streaming_array(name_);
  if (streaming_particles_)
  {
    Ps.turnOnPerParticleSK();
    V_sample = tm.checkout_real<1>(name_, Ps);
    if (!is_active)
      evaluate_sp(Ps);
  }
}

void CoulombPBCAA::deleteParticleQuantities()
{
  if (streaming_particles_)
    delete V_sample;
}
#endif

void CoulombPBCAA::informOfPerParticleListener()
{
  // turnOnParticleSK is written so it can be called again and again.
  Ps.turnOnPerParticleSK();
  OperatorBase::informOfPerParticleListener();
}


CoulombPBCAA::Return_t CoulombPBCAA::evaluate(ParticleSet& P)
{
  if (is_active)
  {
#if !defined(REMOVE_TRACEMANAGER)
    if (streaming_particles_)
      value_ = evaluate_sp(P);
    else
#endif
      value_ = evalLR(P) + evalSR(P) + myConst;
  }
  return value_;
}

void CoulombPBCAA::mw_evaluate(const RefVectorWithLeader<OperatorBase>& o_list,
                               const RefVectorWithLeader<TrialWaveFunction>& wf_list,
                               const RefVectorWithLeader<ParticleSet>& p_list) const
{
  auto& o_leader = o_list.getCastedLeader<CoulombPBCAA>();
  auto& p_leader = p_list.getLeader();
  assert(this == &o_list.getLeader());

  if (!o_leader.is_active)
    return;

  if (use_offload_)
  {
    if (o_leader.streaming_particles_)
      throw std::runtime_error("Streaming particles is not supported when offloading in CoulombPBCAA");

    auto short_range_results = mw_evalSR_offload(o_list, p_list);

    for (int iw = 0; iw < o_list.size(); iw++)
    {
      auto& coulomb_aa  = o_list.getCastedElement<CoulombPBCAA>(iw);
      coulomb_aa.value_ = coulomb_aa.evalLR(p_list[iw]) + short_range_results[iw] + myConst;
    }
  }
  else
    OperatorBase::mw_evaluate(o_list, wf_list, p_list);
}

void CoulombPBCAA::mw_evaluatePerParticle(const RefVectorWithLeader<OperatorBase>& o_list,
                                          const RefVectorWithLeader<TrialWaveFunction>& wf_list,
                                          const RefVectorWithLeader<ParticleSet>& p_list,
                                          const std::vector<ListenerVector<RealType>>& listeners,
                                          const std::vector<ListenerVector<RealType>>& listeners_ions) const
{
  auto& o_leader = o_list.getCastedLeader<CoulombPBCAA>();
  auto& p_leader = p_list.getLeader();
  assert(this == &o_list.getLeader());

  if (!o_leader.is_active)
    return;

  auto num_centers = p_leader.getTotalNum();
  auto name(o_leader.getName());
  Vector<RealType>& v_sample = o_leader.mw_res_handle_.getResource().v_sample;
  const auto& pp_consts      = o_leader.mw_res_handle_.getResource().pp_consts;
  auto num_species           = p_leader.getSpeciesSet().getTotalNum();
  v_sample.resize(num_centers);
  // This lambda is mostly about getting a handle on what is being touched by the per particle evaluation.
  auto evaluate_walker = [num_species, num_centers, name, &v_sample,
                          &pp_consts](const int walker_index, const CoulombPBCAA& cpbcaa, const ParticleSet& pset,
                                      const std::vector<ListenerVector<RealType>>& listeners) -> RealType {
    mRealType Vsr = 0.0;
    mRealType Vlr = 0.0;
    mRealType Vc  = cpbcaa.myConst;
    std::fill(v_sample.begin(), v_sample.end(), 0.0);
    {
      //SR
      const auto& d_aa(pset.getDistTableAA(cpbcaa.d_aa_ID));
      RealType z;
      for (int ipart = 1; ipart < num_centers; ipart++)
      {
        z                = .5 * cpbcaa.Zat[ipart];
        const auto& dist = d_aa.getDistRow(ipart);
        for (int jpart = 0; jpart < ipart; ++jpart)
        {
          RealType pairpot = z * cpbcaa.Zat[jpart] * cpbcaa.rVs->splint(dist[jpart]) / dist[jpart];
          v_sample[ipart] += pairpot;
          v_sample[jpart] += pairpot;
          Vsr += pairpot;
        }
      }
      Vsr *= 2.0;
    }
    {
      //LR
      const StructFact& PtclRhoK(pset.getSK());
      if (PtclRhoK.SuperCellEnum == SUPERCELL_SLAB)
      {
        APP_ABORT("CoulombPBCAA::evaluate_sp single particle traces have not been implemented for slab geometry");
      }
      else
      {
        assert(PtclRhoK.isStorePerParticle()); // ensure this so we know eikr_r has been allocated
        //jtk mark: needs optimizations
        RealType v1; //single particle energy
        RealType z;
        for (int i = 0; i < num_centers; i++)
        {
          z  = .5 * cpbcaa.Zat[i];
          v1 = 0.0;
          for (int s = 0; s < num_species; ++s)
            v1 += z * cpbcaa.Zspec[s] *
                cpbcaa.AA->evaluate(pset.getSimulationCell().getKLists().kshell, PtclRhoK.rhok_r[s], PtclRhoK.rhok_i[s],
                                    PtclRhoK.eikr_r[i], PtclRhoK.eikr_i[i]);
          v_sample[i] += v1;
          Vlr += v1;
        }
      }
    }
    for (int i = 0; i < v_sample.size(); ++i)
      v_sample[i] += pp_consts[i];
    RealType value = Vsr + Vlr + Vc;

    for (const ListenerVector<RealType>& listener : listeners)
      listener.report(walker_index, name, v_sample);
    return value;
  };

  for (int iw = 0; iw < o_list.size(); iw++)
  {
    auto& coulomb_aa  = o_list.getCastedElement<CoulombPBCAA>(iw);
    coulomb_aa.value_ = evaluate_walker(iw, coulomb_aa, p_list[iw], listeners);
  }
}

CoulombPBCAA::Return_t CoulombPBCAA::evaluateWithIonDerivs(ParticleSet& P,
                                                           ParticleSet& ions,
                                                           TrialWaveFunction& psi,
                                                           ParticleSet::ParticlePos& hf_terms,
                                                           ParticleSet::ParticlePos& pulay_terms)
{
  if (ComputeForces and !is_active)
    hf_terms -= forces_;
  //No pulay term.
  return value_;
}

#if !defined(REMOVE_TRACEMANAGER)
CoulombPBCAA::Return_t CoulombPBCAA::evaluate_sp(ParticleSet& P)
{
  mRealType Vsr              = 0.0;
  mRealType Vlr              = 0.0;
  mRealType& Vc              = myConst;
  Array<RealType, 1>& V_samp = *V_sample;
  V_samp                     = 0.0;
  {
    //SR
    const auto& d_aa(P.getDistTableAA(d_aa_ID));
    RealType z;
    for (int ipart = 1; ipart < NumCenters; ipart++)
    {
      z                = .5 * Zat[ipart];
      const auto& dist = d_aa.getDistRow(ipart);
      for (int jpart = 0; jpart < ipart; ++jpart)
      {
        RealType pairpot = z * Zat[jpart] * rVs->splint(dist[jpart]) / dist[jpart];
        V_samp(ipart) += pairpot;
        V_samp(jpart) += pairpot;
        Vsr += pairpot;
      }
    }
    Vsr *= 2.0;
  }
  {
    //LR
    const StructFact& PtclRhoK(P.getSK());
    if (PtclRhoK.SuperCellEnum == SUPERCELL_SLAB)
    {
      APP_ABORT("CoulombPBCAA::evaluate_sp single particle traces have not been implemented for slab geometry");
    }
    else
    {
      assert(PtclRhoK.isStorePerParticle()); // ensure this so we know eikr_r has been allocated
      //jtk mark: needs optimizations
      RealType v1; //single particle energy
      RealType z;
      for (int i = 0; i < NumCenters; i++)
      {
        z  = .5 * Zat[i];
        v1 = 0.0;
        for (int s = 0; s < NumSpecies; ++s)
          v1 += z * Zspec[s] *
              AA->evaluate(P.getSimulationCell().getKLists().kshell, PtclRhoK.rhok_r[s], PtclRhoK.rhok_i[s],
                           PtclRhoK.eikr_r[i], PtclRhoK.eikr_i[i]);
        V_samp(i) += v1;
        Vlr += v1;
      }
    }
  }
  for (int i = 0; i < V_samp.size(); ++i)
    V_samp(i) += V_const(i);
  value_ = Vsr + Vlr + Vc;
#if defined(TRACE_CHECK)
  RealType Vlrnow = evalLR(P);
  RealType Vsrnow = evalSR(P);
  RealType Vcnow  = myConst;
  RealType Vnow   = Vlrnow + Vsrnow + Vcnow;
  RealType Vsum   = V_samp.sum();
  RealType Vcsum  = V_const.sum();
  if (std::abs(Vsum - Vnow) > TraceManager::trace_tol)
  {
    app_log() << "accumtest: CoulombPBCAA::evaluate()" << std::endl;
    app_log() << "accumtest:   tot:" << Vnow << std::endl;
    app_log() << "accumtest:   sum:" << Vsum << std::endl;
    APP_ABORT("Trace check failed");
  }
  if (std::abs(Vcsum - Vcnow) > TraceManager::trace_tol)
  {
    app_log() << "accumtest: CoulombPBCAA::evalConsts()" << std::endl;
    app_log() << "accumtest:   tot:" << Vcnow << std::endl;
    app_log() << "accumtest:   sum:" << Vcsum << std::endl;
    APP_ABORT("Trace check failed");
  }
#endif
  return value_;
}
#endif

void CoulombPBCAA::initBreakup(ParticleSet& P)
{
  //SpeciesSet& tspecies(PtclRef->getSpeciesSet());
  SpeciesSet& tspecies(P.getSpeciesSet());
  //Things that don't change with lattice are done here instead of InitBreakup()
  ChargeAttribIndx = tspecies.addAttribute("charge");
  NumCenters       = P.getTotalNum();
  NumSpecies       = tspecies.TotalNum;

#if !defined(REMOVE_TRACEMANAGER)
  V_const.resize(NumCenters);
#endif

  Zspec.resize(NumSpecies);
  NofSpecies.resize(NumSpecies);
  for (int spec = 0; spec < NumSpecies; spec++)
  {
    Zspec[spec]      = tspecies(ChargeAttribIndx, spec);
    NofSpecies[spec] = P.groupsize(spec);
  }

  SpeciesID.resize(NumCenters);
  Zat.resize(NumCenters);
  Zat_offload = std::make_shared<Vector<RealType, OffloadPinnedAllocator<RealType>>>(NumCenters);
  auto& Zat_ref(*Zat_offload);
  for (int iat = 0; iat < NumCenters; iat++)
  {
    SpeciesID[iat] = P.GroupID[iat];
    Zat[iat]       = Zspec[P.GroupID[iat]];
    Zat_ref[iat]   = Zat[iat];
  }
  Zat_ref.updateTo();

  AA = LRCoulombSingleton::getHandler(P);
  //AA->initBreakup(*PtclRef);
  myConst = evalConsts();
  myRcut  = AA->get_rc(); //Basis.get_rc();

  if (rVs == nullptr)
    rVs = LRCoulombSingleton::createSpline4RbyVs(AA.get(), myRcut);

  rVs_offload = std::make_shared<const OffloadSpline>(*rVs);

  if (ComputeForces)
  {
    dAA = LRCoulombSingleton::getDerivHandler(P);
    if (rVsforce == nullptr)
    {
      rVsforce = LRCoulombSingleton::createSpline4RbyVs(dAA.get(), myRcut);
    }
  }
}


CoulombPBCAA::Return_t CoulombPBCAA::evalLRwithForces(ParticleSet& P)
{
  //  const StructFact& PtclRhoK(P.getSK());
  std::vector<TinyVector<RealType, DIM>> grad(P.getTotalNum());
  for (int spec2 = 0; spec2 < NumSpecies; spec2++)
  {
    RealType Z2 = Zspec[spec2];
    for (int iat = 0; iat < grad.size(); iat++)
      grad[iat] = TinyVector<RealType, DIM>(0.0);
    //AA->evaluateGrad(P, P, spec2, Zat, grad);
    dAA->evaluateGrad(P, P, spec2, Zat, grad);
    for (int iat = 0; iat < grad.size(); iat++)
      forces_[iat] += Z2 * grad[iat];
  } //spec2
  return evalLR(P);
}


CoulombPBCAA::Return_t CoulombPBCAA::evalSRwithForces(ParticleSet& P)
{
  const auto& d_aa(P.getDistTableAA(d_aa_ID));
  mRealType SR = 0.0;
  for (size_t ipart = 1; ipart < (NumCenters / 2 + 1); ipart++)
  {
    mRealType esum   = 0.0;
    const auto& dist = d_aa.getDistRow(ipart);
    const auto& dr   = d_aa.getDisplRow(ipart);
    for (size_t j = 0; j < ipart; ++j)
    {
      RealType V, rV, d_rV_dr, d2_rV_dr2;
      RealType rinv = 1.0 / dist[j];
      rV            = rVsforce->splint(dist[j], d_rV_dr, d2_rV_dr2);
      V             = rV * rinv;
      esum += Zat[j] * rVs->splint(dist[j]) * rinv;

      PosType grad = Zat[j] * Zat[ipart] * (d_rV_dr - V) * rinv * rinv * dr[j];
      forces_[ipart] += grad;
      forces_[j] -= grad;
    }
    SR += Zat[ipart] * esum;

    const size_t ipart_reverse = NumCenters - ipart;
    if (ipart == ipart_reverse)
      continue;

    esum              = 0.0;
    const auto& dist2 = d_aa.getDistRow(ipart_reverse);
    const auto& dr2   = d_aa.getDisplRow(ipart_reverse);
    for (size_t j = 0; j < ipart_reverse; ++j)
    {
      RealType V, rV, d_rV_dr, d2_rV_dr2;
      RealType rinv = 1.0 / dist2[j];
      rV            = rVsforce->splint(dist2[j], d_rV_dr, d2_rV_dr2);
      V             = rV * rinv;
      esum += Zat[j] * rVs->splint(dist2[j]) * rinv;

      PosType grad = Zat[j] * Zat[ipart_reverse] * (d_rV_dr - V) * rinv * rinv * dr2[j];
      forces_[ipart_reverse] += grad;
      forces_[j] -= grad;
    }
    SR += Zat[ipart_reverse] * esum;
  }
  return SR;
}


/** evaluate the constant term that does not depend on the position
 *
 * \htmlonly
 * <ul>
 * <li> self-energy: \f$ -\frac{1}{2}\sum_{i} v_l(r=0) q_i^2 = -\frac{1}{2}v_l(r=0) \sum_{alpha} N^{\alpha} q^{\alpha}^2\f$
 * <li> background term \f$ V_{bg} = -\frac{1}{2}\sum_{\alpha}\sum_{\beta} N^{\alpha}q^{\alpha}N^{\beta}q^{\beta} v_s(k=0)\f$
 * </ul>
 * \endhtmlonly
 * CoulombPBCABTemp contributes additional background term which completes the background term
 * note this calculates the per particle consts even if trace manager is removed.
 */
CoulombPBCAA::Return_t CoulombPBCAA::evalConsts(bool report)
{
  mRealType Consts = 0.0; // constant term
  mRealType v1;           //single particle energy
  mRealType vl_r0 = AA->evaluateLR_r0();
  mRealType vs_k0 = AA->evaluateSR_k0();

  if (quasi2d) // background term has z dependence
  {            // just evaluate the Madelung term
    for (int ispec = 1; ispec < NumSpecies; ispec++)
      if (Zspec[ispec] != Zspec[0])
        throw std::runtime_error("quasi2d assumes same charge");
    if (report)
    {
      app_log() << "    vlr(r->0) = " << vl_r0 << std::endl;
      app_log() << "   1/V vsr_k0 = " << vs_k0 << std::endl;
    }
    // make sure we can ignore the short-range Madelung sum
    mRealType Rws           = Ps.getLattice().WignerSeitzRadius;
    mRealType rvsr_at_image = Rws * AA->evaluate(Rws, 1.0 / Rws);
    if (rvsr_at_image > 1e-6)
    {
      std::ostringstream msg;
      msg << std::setprecision(14);
      msg << "Ewald alpha = " << rvsr_at_image << " is too small" << std::endl;
      msg << "Short-range potential r*vsr(r) = " << rvsr_at_image << " at image radius r=" << Rws << std::endl;
      throw std::runtime_error(msg.str());
    }
    // perform long-range Madelung sum
    const StructFact& PtclRhoK(Ps.getSK());
    v1 = AA->evaluate_slab(0, Ps.getSimulationCell().getKLists().kshell, PtclRhoK.eikr_r[0], PtclRhoK.eikr_i[0],
                           PtclRhoK.eikr_r[0], PtclRhoK.eikr_i[0]);
    if (report)
      app_log() << "   LR Madelung = " << v1 << std::endl;
    madelung_constant_ = 0.5 * (v1 - vl_r0);
    Consts             = NumCenters * madelung_constant_;
  }
  else // group background term together with Madelung vsr_k0 part
  {
#if !defined(REMOVE_TRACEMANAGER)
    V_const = 0.0;
#endif
    for (int ipart = 0; ipart < NumCenters; ipart++)
    {
      v1 = -.5 * Zat[ipart] * Zat[ipart] * vl_r0;
#if !defined(REMOVE_TRACEMANAGER)
      V_const(ipart) += v1;
#endif
      Consts += v1;
    }
    if (report)
      app_log() << "   PBCAA self-interaction term " << Consts << std::endl;
    //Compute Madelung constant
    madelung_constant_ = 0.0;
    for (int i = 0; i < AA->Fk.size(); i++)
      madelung_constant_ += AA->Fk[i];
    madelung_constant_ = 0.5 * (madelung_constant_ - vl_r0 - vs_k0);
    for (int ipart = 0; ipart < NumCenters; ipart++)
    {
      v1 = 0.0;
      for (int spec = 0; spec < NumSpecies; spec++)
        v1 += NofSpecies[spec] * Zspec[spec];
      v1 *= -.5 * Zat[ipart] * vs_k0;
#if !defined(REMOVE_TRACEMANAGER)
      V_const(ipart) += v1;
#endif
      Consts += v1;
    }
    if (report)
      app_log() << "   PBCAA total constant " << Consts << std::endl;
  }
  return Consts;
}


CoulombPBCAA::Return_t CoulombPBCAA::evalSR(ParticleSet& P)
{
  ScopedTimer local_timer(evalSR_timer_);
  const auto& d_aa(P.getDistTableAA(d_aa_ID));
  mRealType SR = 0.0;
#pragma omp parallel for reduction(+ : SR)
  for (size_t ipart = 1; ipart < (NumCenters / 2 + 1); ipart++)
  {
    mRealType esum   = 0.0;
    const auto& dist = d_aa.getDistRow(ipart);
    for (size_t j = 0; j < ipart; ++j)
      esum += Zat[j] * rVs->splint(dist[j]) / dist[j];
    SR += Zat[ipart] * esum;

    const size_t ipart_reverse = NumCenters - ipart;
    if (ipart == ipart_reverse)
      continue;

    esum              = 0.0;
    const auto& dist2 = d_aa.getDistRow(ipart_reverse);
    for (size_t j = 0; j < ipart_reverse; ++j)
      esum += Zat[j] * rVs->splint(dist2[j]) / dist2[j];
    SR += Zat[ipart_reverse] * esum;
  }
  return SR;
}

std::vector<CoulombPBCAA::Return_t> CoulombPBCAA::mw_evalSR_offload(const RefVectorWithLeader<OperatorBase>& o_list,
                                                                    const RefVectorWithLeader<ParticleSet>& p_list)
{
  const size_t nw  = o_list.size();
  auto& p_leader   = p_list.getLeader();
  auto& caa_leader = o_list.getCastedLeader<CoulombPBCAA>();
  ScopedTimer local_timer(caa_leader.evalSR_timer_);

  RefVectorWithLeader<DistanceTable> dt_list(p_leader.getDistTable(caa_leader.d_aa_ID));
  dt_list.reserve(p_list.size());
  for (ParticleSet& p : p_list)
    dt_list.push_back(p.getDistTable(caa_leader.d_aa_ID));

  auto& dtaa_leader = dynamic_cast<DistanceTableAA&>(p_leader.getDistTable(caa_leader.d_aa_ID));

  const size_t chunk_size = dtaa_leader.get_num_particls_stored();
  if (chunk_size == 0)
    throw std::runtime_error("bug dtaa_leader.get_num_particls_stored() == 0");

  auto& values_offload        = caa_leader.mw_res_handle_.getResource().values_offload;
  const size_t total_num      = p_leader.getTotalNum();
  const size_t total_num_half = (total_num + 1) / 2;
  const size_t num_padded     = getAlignedSize<RealType>(total_num);
  const size_t num_chunks     = (total_num_half + chunk_size - 1) / chunk_size;

  const auto m_Y         = caa_leader.rVs_offload->get_m_Y().data();
  const auto m_Y2        = caa_leader.rVs_offload->get_m_Y2().data();
  const auto first_deriv = caa_leader.rVs_offload->get_first_deriv();
  const auto const_value = caa_leader.rVs_offload->get_const_value();
  const auto r_min       = caa_leader.rVs_offload->get_r_min();
  const auto r_max       = caa_leader.rVs_offload->get_r_max();
  const auto X           = caa_leader.rVs_offload->get_X().data();
  const auto delta_inv   = caa_leader.rVs_offload->get_delta_inv();
  const auto Zat         = caa_leader.Zat_offload->data();

  {
    values_offload.resize(nw);
    std::fill_n(values_offload.data(), nw, 0);
    auto value_ptr = values_offload.data();
    values_offload.updateTo();
    for (size_t ichunk = 0; ichunk < num_chunks; ichunk++)
    {
      const size_t first           = ichunk * chunk_size;
      const size_t last            = std::min(first + chunk_size, total_num_half);
      const size_t this_chunk_size = last - first;

      auto* mw_dist = dtaa_leader.mw_evalDistsInRange(dt_list, p_list, first, last);

      ScopedTimer offload_scope(caa_leader.offload_timer_);

      PRAGMA_OFFLOAD("omp target teams distribute num_teams(nw)")
      for (size_t iw = 0; iw < nw; iw++)
      {
        mRealType SR = 0.0;
        PRAGMA_OFFLOAD("omp parallel for reduction(+ : SR)")
        for (size_t jcol = 0; jcol < total_num; jcol++)
          for (size_t irow = first; irow < last; irow++)
          {
            const RealType dist = mw_dist[num_padded * (irow - first + iw * this_chunk_size) + jcol];
            if (irow == jcol || (irow * 2 + 1 == total_num && jcol > irow))
              continue;

            const size_t i = irow > jcol ? irow : total_num - 1 - irow;
            const size_t j = irow > jcol ? jcol : total_num - 1 - jcol;

            SR += Zat[i] * Zat[j] *
                OffloadSpline::splint(r_min, r_max, X, delta_inv, m_Y, m_Y2, first_deriv, const_value, dist) / dist;
          }
        value_ptr[iw] += SR;
      }
    }

    values_offload.updateFrom();
  }
  std::vector<Return_t> values(nw);
  for (int iw = 0; iw < nw; iw++)
    values[iw] = values_offload[iw];
  return values;
}

CoulombPBCAA::Return_t CoulombPBCAA::evalLR(ParticleSet& P)
{
  ScopedTimer local_timer(evalLR_timer_);
  mRealType res = 0.0;
  const StructFact& PtclRhoK(P.getSK());
  if (quasi2d)
  {
    const auto& d_aa(P.getDistTableAA(d_aa_ID));
    // need 1/2 \sum_{i,j} v_E(r_i - r_j)
    //distance table handles jat<iat
    for (int iat = 1; iat < NumCenters; ++iat)
    {
      mRealType u        = 0;
      const int slab_dir = OHMMS_DIM - 1;
      const auto& dr     = d_aa.getDisplRow(iat);
      for (int jat = 0; jat < iat; ++jat)
      {
        const RealType z = std::abs(dr[jat][slab_dir]);
        u += Zat[jat] *
            AA->evaluate_slab(z, P.getSimulationCell().getKLists().kshell, PtclRhoK.eikr_r[iat], PtclRhoK.eikr_i[iat],
                              PtclRhoK.eikr_r[jat], PtclRhoK.eikr_i[jat]);
      }
      res += Zat[iat] * u;
    }
  }
  else
  {
    for (int spec1 = 0; spec1 < NumSpecies; spec1++)
    {
      mRealType Z1 = Zspec[spec1];
      for (int spec2 = spec1; spec2 < NumSpecies; spec2++)
      {
        mRealType temp = AA->evaluate(P.getSimulationCell().getKLists().kshell, PtclRhoK.rhok_r[spec1],
                                      PtclRhoK.rhok_i[spec1], PtclRhoK.rhok_r[spec2], PtclRhoK.rhok_i[spec2]);
        if (spec2 == spec1)
          temp *= 0.5;
        res += Z1 * Zspec[spec2] * temp;
      } //spec2
    }   //spec1
  }
  return res;
}

void CoulombPBCAA::evalPerParticleConsts(Vector<RealType>& pp_consts) const
{
  mRealType v1; //single particle energy
  mRealType vl_r0 = AA->evaluateLR_r0();
  mRealType vs_k0 = AA->evaluateSR_k0();

  if (quasi2d)
    throw std::runtime_error("Batched per particle eval is not supported for quasi2d");
  else
  {
    pp_consts.resize(NumCenters, 0.0);
    for (int ipart = 0; ipart < NumCenters; ipart++)
    {
      v1 = -.5 * Zat[ipart] * Zat[ipart] * vl_r0;
      pp_consts[ipart] += v1;
    }
    for (int ipart = 0; ipart < NumCenters; ipart++)
    {
      v1 = 0.0;
      for (int spec = 0; spec < NumSpecies; spec++)
        v1 += NofSpecies[spec] * Zspec[spec];
      v1 *= -.5 * Zat[ipart] * vs_k0;
      pp_consts[ipart] += v1;
    }
  }
}

void CoulombPBCAA::createResource(ResourceCollection& collection) const
{
  auto new_res = std::make_unique<CoulombPBCAAMultiWalkerResource>();
  if (hasListener())
    evalPerParticleConsts(new_res->pp_consts);
  auto resource_index = collection.addResource(std::move(new_res));
}

void CoulombPBCAA::acquireResource(ResourceCollection& collection,
                                   const RefVectorWithLeader<OperatorBase>& o_list) const
{
  auto& o_leader          = o_list.getCastedLeader<CoulombPBCAA>();
  o_leader.mw_res_handle_ = collection.lendResource<CoulombPBCAAMultiWalkerResource>();
}

void CoulombPBCAA::releaseResource(ResourceCollection& collection,
                                   const RefVectorWithLeader<OperatorBase>& o_list) const
{
  auto& o_leader = o_list.getCastedLeader<CoulombPBCAA>();
  collection.takebackResource(o_leader.mw_res_handle_);
}

std::unique_ptr<OperatorBase> CoulombPBCAA::makeClone(ParticleSet& qp, TrialWaveFunction& psi)
{
  return std::make_unique<CoulombPBCAA>(*this);
}
} // namespace qmcplusplus
