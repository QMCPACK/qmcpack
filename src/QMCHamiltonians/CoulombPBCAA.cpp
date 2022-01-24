//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "EwaldRef.h"
#include "CoulombPBCAA.h"
#include "Particle/DistanceTable.h"
#include "Utilities/ProgressReportEngine.h"
#include <numeric>

namespace qmcplusplus
{
CoulombPBCAA::CoulombPBCAA(ParticleSet& ref, bool active, bool computeForces)
    : ForceBase(ref, ref),
      is_active(active),
      FirstTime(true),
      myConst(0.0),
      ComputeForces(computeForces),
      Ps(ref),
      d_aa_ID(ref.addTable(ref)),
      evalLR_timer_(*timer_manager.createTimer("CoulombPBCAA::LongRange", timer_level_fine)),
      evalSR_timer_(*timer_manager.createTimer("CoulombPBCAA::ShortRange", timer_level_fine))

{
  ReportEngine PRE("CoulombPBCAA", "CoulombPBCAA");
  setEnergyDomain(POTENTIAL);
  twoBodyQuantumDomain(ref);
  PtclRefName = ref.getDistTable(d_aa_ID).getName();
  initBreakup(ref);

  if (ComputeForces)
  {
    ref.turnOnPerParticleSK();
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

      throw std::runtime_error(msg.str());
    }
    else
    {
      app_log() << "  Check passed." << std::endl;
    }
  }
  prefix = "F_AA";
  app_log() << "  Maximum K shell " << AA->MaxKshell << std::endl;
  app_log() << "  Number of k vectors " << AA->Fk.size() << std::endl;
  app_log() << "  Fixed Coulomb potential for " << ref.getName();
  app_log() << "\n    e-e Madelung Const. =" << MC0 << "\n    Vtot     =" << value_ << std::endl;
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
    forces = 0.0;
    eS     = evalSRwithForces(s);
    eL     = evalLRwithForces(s);
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

CoulombPBCAA::Return_t CoulombPBCAA::evaluateWithIonDerivs(ParticleSet& P,
                                                           ParticleSet& ions,
                                                           TrialWaveFunction& psi,
                                                           ParticleSet::ParticlePos& hf_terms,
                                                           ParticleSet::ParticlePos& pulay_terms)
{
  if (ComputeForces and !is_active)
    hf_terms -= forces;
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
      //jtk mark: needs optimizations for USE_REAL_STRUCT_FACTOR
      RealType v1; //single particle energy
      RealType z;
      for (int i = 0; i < NumCenters; i++)
      {
        z  = .5 * Zat[i];
        v1 = 0.0;
        for (int s = 0; s < NumSpecies; ++s)
        {
#if defined(USE_REAL_STRUCT_FACTOR)
          v1 += z * Zspec[s] *
              AA->evaluate(PtclRhoK.getKLists().kshell, PtclRhoK.rhok_r[s], PtclRhoK.rhok_i[s], PtclRhoK.eikr_r[i],
                           PtclRhoK.eikr_i[i]);
#else
          v1 += z * Zspec[s] * AA->evaluate(PtclRhoK.getKLists().kshell, PtclRhoK.rhok[s], PtclRhoK.eikr[i]);
#endif
        }
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
  MemberAttribIndx = tspecies.addAttribute("membersize");
  NumCenters       = P.getTotalNum();
  NumSpecies       = tspecies.TotalNum;

#if !defined(REMOVE_TRACEMANAGER)
  V_const.resize(NumCenters);
#endif

  Zat.resize(NumCenters);
  Zspec.resize(NumSpecies);
  NofSpecies.resize(NumSpecies);
  for (int spec = 0; spec < NumSpecies; spec++)
  {
    Zspec[spec]      = tspecies(ChargeAttribIndx, spec);
    NofSpecies[spec] = static_cast<int>(tspecies(MemberAttribIndx, spec));
  }
  SpeciesID.resize(NumCenters);
  for (int iat = 0; iat < NumCenters; iat++)
  {
    SpeciesID[iat] = P.GroupID[iat];
    Zat[iat]       = Zspec[P.GroupID[iat]];
  }
  AA = LRCoulombSingleton::getHandler(P);
  //AA->initBreakup(*PtclRef);
  myConst = evalConsts();
  myRcut  = AA->get_rc(); //Basis.get_rc();

  if (rVs == nullptr)
    rVs = LRCoulombSingleton::createSpline4RbyVs(AA.get(), myRcut);

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
      forces[iat] += Z2 * grad[iat];
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
      forces[ipart] += grad;
      forces[j] -= grad;
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
      forces[ipart_reverse] += grad;
      forces[j] -= grad;
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
 */
CoulombPBCAA::Return_t CoulombPBCAA::evalConsts(bool report)
{
  mRealType Consts = 0.0; // constant term
  mRealType v1;           //single particle energy
#if !defined(REMOVE_TRACEMANAGER)
  V_const = 0.0;
#endif
  //v_l(r=0) including correction due to the non-periodic direction
  mRealType vl_r0 = AA->evaluateLR_r0();
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
  //Compute Madelung constant: this is not correct for general cases
  MC0 = 0.0;
  for (int i = 0; i < AA->Fk.size(); i++)
    MC0 += AA->Fk[i];
  MC0 = 0.5 * (MC0 - vl_r0);
  //Neutraling background term
  mRealType vs_k0 = AA->evaluateSR_k0(); //v_s(k=0)
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
  //app_log() << "   MC0 of PBCAA " << MC0 << std::endl;
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

CoulombPBCAA::Return_t CoulombPBCAA::evalLR(ParticleSet& P)
{
  ScopedTimer local_timer(evalLR_timer_);
  mRealType res = 0.0;
  const StructFact& PtclRhoK(P.getSK());
  if (PtclRhoK.SuperCellEnum == SUPERCELL_SLAB)
  {
    const auto& d_aa(P.getDistTableAA(d_aa_ID));
    //distance table handles jat<iat
    for (int iat = 1; iat < NumCenters; ++iat)
    {
      mRealType u = 0;
#if !defined(USE_REAL_STRUCT_FACTOR)
      const int slab_dir              = OHMMS_DIM - 1;
      const RealType* restrict d_slab = d_aa.Displacements[iat].data(slab_dir);
      for (int jat = 0; jat < iat; ++jat)
        u += Zat[jat] *
            AA->evaluate_slab(-d_slab[jat], //JK: Could be wrong. Check the SIGN
                              PtclRhoK.getKLists().kshell, PtclRhoK.eikr[iat], PtclRhoK.eikr[jat]);
#endif
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
#if defined(USE_REAL_STRUCT_FACTOR)
        mRealType temp = AA->evaluate(PtclRhoK.getKLists().kshell, PtclRhoK.rhok_r[spec1], PtclRhoK.rhok_i[spec1],
                                      PtclRhoK.rhok_r[spec2], PtclRhoK.rhok_i[spec2]);
#else
        mRealType temp = AA->evaluate(PtclRhoK.getKLists().kshell, PtclRhoK.rhok[spec1], PtclRhoK.rhok[spec2]);
#endif
        if (spec2 == spec1)
          temp *= 0.5;
        res += Z1 * Zspec[spec2] * temp;
      } //spec2
    }   //spec1
  }
  return res;
}

std::unique_ptr<OperatorBase> CoulombPBCAA::makeClone(ParticleSet& qp, TrialWaveFunction& psi)
{
  return std::make_unique<CoulombPBCAA>(*this);
}
} // namespace qmcplusplus
