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

#include <ResourceCollection.h>
#include "CoulombPBCAB.h"
#include "Particle/DistanceTable.h"
#include "Message/Communicate.h"
#include "Utilities/ProgressReportEngine.h"

namespace qmcplusplus
{
struct CoulombPBCAB::CoulombPBCABMultiWalkerResource : public Resource
{
  CoulombPBCABMultiWalkerResource() : Resource("CoulombPBCAB") {}
  std::unique_ptr<Resource> makeClone() const override
  {
    return std::make_unique<CoulombPBCABMultiWalkerResource>(*this);
  }

  /// a walkers worth of per ion AB potential values
  Vector<RealType> pp_samples_src;
  /// a walkers worth of per electron AB potential values
  Vector<RealType> pp_samples_trg;
  /// constant values for the source particles aka ions aka A
  Vector<RealType> pp_consts_src;
  /// constant values for the target particles aka electrons aka B
  Vector<RealType> pp_consts_trg;
};

CoulombPBCAB::CoulombPBCAB(ParticleSet& ions, ParticleSet& elns, bool computeForces)
    : ForceBase(ions, elns),
      myTableIndex(elns.addTable(ions)),
      myConst(0.0),
      ComputeForces(computeForces),
      MaxGridPoints(10000),
      Peln(elns),
      pset_ions_(ions)
{
  ReportEngine PRE("CoulombPBCAB", "CoulombPBCAB");
  setEnergyDomain(POTENTIAL);
  twoBodyQuantumDomain(ions, elns);
  if (ComputeForces)
    pset_ions_.turnOnPerParticleSK();
  initBreakup(elns);
  prefix_ = "Flocal";
  app_log() << "  Rcut                " << myRcut << std::endl;
  app_log() << "  Maximum K shell     " << AB->MaxKshell << std::endl;
  app_log() << "  Number of k vectors " << AB->Fk.size() << std::endl;
}

std::unique_ptr<OperatorBase> CoulombPBCAB::makeClone(ParticleSet& qp, TrialWaveFunction& psi)
{
  return std::make_unique<CoulombPBCAB>(*this);
}

CoulombPBCAB::~CoulombPBCAB() = default;

void CoulombPBCAB::resetTargetParticleSet(ParticleSet& P)
{
  int tid = P.addTable(pset_ions_);
  if (tid != myTableIndex)
  {
    APP_ABORT("CoulombPBCAB::resetTargetParticleSet found inconsistent table index");
  }
  AB->resetTargetParticleSet(P);
}

void CoulombPBCAB::addObservables(PropertySetType& plist, BufferType& collectables)
{
  my_index_ = plist.add(name_.c_str());
  if (ComputeForces)
    addObservablesF(plist);
}


#if !defined(REMOVE_TRACEMANAGER)
void CoulombPBCAB::contributeParticleQuantities() { request_.contribute_array(name_); }

void CoulombPBCAB::checkoutParticleQuantities(TraceManager& tm)
{
  streaming_particles_ = request_.streaming_array(name_);
  if (streaming_particles_)
  {
    pset_ions_.turnOnPerParticleSK();
    Peln.turnOnPerParticleSK();
    Ve_sample = tm.checkout_real<1>(name_, Peln);
    Vi_sample = tm.checkout_real<1>(name_, pset_ions_);
  }
}

void CoulombPBCAB::informOfPerParticleListener()
{
  // This is written so it can be called again and again.
  pset_ions_.turnOnPerParticleSK();
  Peln.turnOnPerParticleSK();
  OperatorBase::informOfPerParticleListener();
}

void CoulombPBCAB::deleteParticleQuantities()
{
  if (streaming_particles_)
  {
    delete Ve_sample;
    delete Vi_sample;
  }
}
#endif


CoulombPBCAB::Return_t CoulombPBCAB::evaluate(ParticleSet& P)
{
  if (ComputeForces)
  {
    forces_ = 0.0;
    value_  = evalLRwithForces(P) + evalSRwithForces(P) + myConst;
  }
  else
#if !defined(REMOVE_TRACEMANAGER)
      if (streaming_particles_)
    value_ = evaluate_sp(P);
  else
#endif
    value_ = evalLR(P) + evalSR(P) + myConst;
  return value_;
}

CoulombPBCAB::Return_t CoulombPBCAB::evaluateWithIonDerivs(ParticleSet& P,
                                                           ParticleSet& ions,
                                                           TrialWaveFunction& psi,
                                                           ParticleSet::ParticlePos& hf_terms,
                                                           ParticleSet::ParticlePos& pulay_terms)
{
  if (ComputeForces)
  {
    forces_ = 0.0;
    value_  = evalLRwithForces(P) + evalSRwithForces(P) + myConst;
    hf_terms -= forces_;
    //And no Pulay contribution.
  }
  else
    value_ = evalLR(P) + evalSR(P) + myConst;
  return value_;
}

#if !defined(REMOVE_TRACEMANAGER)
CoulombPBCAB::Return_t CoulombPBCAB::evaluate_sp(ParticleSet& P)
{
  RealType Vsr                = 0.0;
  RealType Vlr                = 0.0;
  mRealType& Vc               = myConst;
  Array<RealType, 1>& Ve_samp = *Ve_sample;
  Array<RealType, 1>& Vi_samp = *Vi_sample;
  Ve_samp                     = 0.0;
  Vi_samp                     = 0.0;
  {
    //SR
    const auto& d_ab(P.getDistTableAB(myTableIndex));
    RealType z;
    //Loop over distinct eln-ion pairs
    for (size_t b = 0; b < NptclB; ++b)
    {
      z                = 0.5 * Qat[b];
      const auto& dist = d_ab.getDistRow(b);
      for (size_t a = 0; a < NptclA; ++a)
      {
        Return_t pairpot = z * Zat[a] * Vat[a]->splint(dist[a]) / dist[a];
        Vi_samp(a) += pairpot;
        Ve_samp(b) += pairpot;
        Vsr += pairpot;
      }
    }
    Vsr *= 2.0;
  }
  {
    //LR
    const StructFact& RhoKA(pset_ions_.getSK());
    const StructFact& RhoKB(P.getSK());
    if (RhoKA.SuperCellEnum == SUPERCELL_SLAB)
    {
      APP_ABORT("CoulombPBCAB::evaluate_sp single particle traces have not been implemented for slab geometry");
    }
    else
    {
      assert(RhoKA.isStorePerParticle()); // ensure this so we know eikr_r has been allocated
      assert(RhoKB.isStorePerParticle());
      //jtk mark: needs optimizations. will likely require new function definitions
      RealType v1; //single particle energy
      RealType q;
      for (int i = 0; i < P.getTotalNum(); ++i)
      {
        q  = .5 * Qat[i];
        v1 = 0.0;
        for (int s = 0; s < NumSpeciesA; s++)
          v1 += Zspec[s] * q *
              AB->evaluate(pset_ions_.getSimulationCell().getKLists().kshell, RhoKA.rhok_r[s], RhoKA.rhok_i[s],
                           RhoKB.eikr_r[i], RhoKB.eikr_i[i]);
        Ve_samp(i) += v1;
        Vlr += v1;
      }
      for (int i = 0; i < pset_ions_.getTotalNum(); ++i)
      {
        q  = .5 * Zat[i];
        v1 = 0.0;
        for (int s = 0; s < NumSpeciesB; s++)
          v1 += Qspec[s] * q *
              AB->evaluate(P.getSimulationCell().getKLists().kshell, RhoKB.rhok_r[s], RhoKB.rhok_i[s], RhoKA.eikr_r[i],
                           RhoKA.eikr_i[i]);
        Vi_samp(i) += v1;
        Vlr += v1;
      }
    }
  }
  for (int i = 0; i < Ve_samp.size(); ++i)
    Ve_samp(i) += Ve_const(i);
  for (int i = 0; i < Vi_samp.size(); ++i)
    Vi_samp(i) += Vi_const(i);
  value_ = Vsr + Vlr + Vc;
#if defined(TRACE_CHECK)
  RealType Vlrnow = evalLR(P);
  RealType Vsrnow = evalSR(P);
  RealType Vcnow  = myConst;
  RealType Vnow   = Vlrnow + Vsrnow + Vcnow;
  RealType Vesum  = Ve_samp.sum();
  RealType Vecsum = Ve_const.sum();
  RealType Visum  = Vi_samp.sum();
  RealType Vicsum = Vi_const.sum();
  RealType Vsum   = Vesum + Visum;
  RealType Vcsum  = Vecsum + Vicsum;
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
  if (std::abs(Vesum - Visum) > TraceManager::trace_tol)
  {
    app_log() << "sharetest: CoulombPBCAB::evaluate()" << std::endl;
    app_log() << "sharetest:   e share:" << Vesum << std::endl;
    app_log() << "sharetest:   i share:" << Visum << std::endl;
  }
  if (std::abs(Vecsum - Vicsum) > TraceManager::trace_tol)
  {
    app_log() << "sharetest: CoulombPBCAB::evalConsts()" << std::endl;
    app_log() << "sharetest:   e share:" << Vecsum << std::endl;
    app_log() << "sharetest:   i share:" << Vicsum << std::endl;
  }

#endif
  return value_;
}
#endif

void CoulombPBCAB::mw_evaluatePerParticle(const RefVectorWithLeader<OperatorBase>& o_list,
                                          const RefVectorWithLeader<TrialWaveFunction>& wf_list,
                                          const RefVectorWithLeader<ParticleSet>& p_list,
                                          const std::vector<ListenerVector<RealType>>& listeners,
                                          const std::vector<ListenerVector<RealType>>& ion_listeners) const
{
  auto& o_leader = o_list.getCastedLeader<CoulombPBCAB>();
  auto& p_leader = p_list.getLeader();
  assert(this == &o_list.getLeader());

  auto num_centers            = p_leader.getTotalNum();
  auto& name                  = o_leader.name_;
  auto& mw_res                = o_leader.mw_res_handle_.getResource();
  Vector<RealType>& ve_sample = mw_res.pp_samples_trg;
  Vector<RealType>& vi_sample = mw_res.pp_samples_src;
  const auto& ve_consts       = mw_res.pp_consts_trg;
  const auto& vi_consts       = mw_res.pp_consts_src;
  auto num_species            = p_leader.getSpeciesSet().getTotalNum();
  ve_sample.resize(NptclB);
  vi_sample.resize(NptclA);

  // This lambda is mostly about getting a handle on what is being touched by the per particle evaluation.
  auto evaluate_walker = [name, &ve_sample, &vi_sample, &ve_consts,
                          &vi_consts](const int walker_index, const CoulombPBCAB& cpbcab, const ParticleSet& pset,
                                      const std::vector<ListenerVector<RealType>>& listeners,
                                      const std::vector<ListenerVector<RealType>>& ion_listeners) -> RealType {
    RealType Vsr = 0.0;
    RealType Vlr = 0.0;
    RealType Vc  = cpbcab.myConst;
    std::fill(ve_sample.begin(), ve_sample.end(), 0.0);
    std::fill(vi_sample.begin(), vi_sample.end(), 0.0);
    auto& pset_source = cpbcab.getSourcePSet();
    {
      //SR
      const auto& d_ab(pset.getDistTableAB(cpbcab.myTableIndex));
      RealType z;
      //Loop over distinct eln-ion pairs
      for (size_t b = 0; b < pset.getTotalNum(); ++b)
      {
        z                = 0.5 * cpbcab.Qat[b];
        const auto& dist = d_ab.getDistRow(b);
        for (size_t a = 0; a < pset_source.getTotalNum(); ++a)
        {
          Return_t pairpot = z * cpbcab.Zat[a] * cpbcab.Vat[a]->splint(dist[a]) / dist[a];
          vi_sample[a] += pairpot;
          ve_sample[b] += pairpot;
          Vsr += pairpot;
        }
      }
      Vsr *= 2.0;
    }
    {
      //LR
      // pset_ions_ is a smelly reference to the ion particle set.
      const StructFact& RhoKA(pset_source.getSK());
      const StructFact& RhoKB(pset.getSK());
      if (RhoKA.SuperCellEnum == SUPERCELL_SLAB)
      {
        APP_ABORT("CoulombPBCAB::evaluate_sp single particle traces have not been implemented for slab geometry");
      }
      else
      {
        assert(RhoKA.isStorePerParticle()); // ensure this so we know eikr_r has been allocated
        assert(RhoKB.isStorePerParticle());
        //jtk mark: needs optimizations. will likely require new function definitions
        RealType v1; //single particle energy
        RealType q;
        int num_species_source = pset_source.getSpeciesSet().size();
        for (int i = 0; i < pset.getTotalNum(); ++i)
        {
          q  = .5 * cpbcab.Qat[i];
          v1 = 0.0;
          for (int s = 0; s < num_species_source; s++)
            v1 += cpbcab.Zspec[s] * q *
                cpbcab.AB->evaluate(pset_source.getSimulationCell().getKLists().kshell, RhoKA.rhok_r[s],
                                    RhoKA.rhok_i[s], RhoKB.eikr_r[i], RhoKB.eikr_i[i]);
          ve_sample[i] += v1;
          Vlr += v1;
        }
        int num_species_target = pset.getSpeciesSet().size();
        for (int i = 0; i < pset_source.getTotalNum(); ++i)
        {
          q  = .5 * cpbcab.Zat[i];
          v1 = 0.0;
          for (int s = 0; s < num_species_target; s++)
            v1 += cpbcab.Qspec[s] * q *
                cpbcab.AB->evaluate(pset.getSimulationCell().getKLists().kshell, RhoKB.rhok_r[s], RhoKB.rhok_i[s],
                                    RhoKA.eikr_r[i], RhoKA.eikr_i[i]);
          vi_sample[i] += v1;
          Vlr += v1;
        }
      }
    }
    for (int i = 0; i < ve_sample.size(); ++i)
      ve_sample[i] += ve_consts[i];
    for (int i = 0; i < vi_sample.size(); ++i)
      vi_sample[i] += vi_consts[i];
    RealType value = Vsr + Vlr + Vc;
    for (const ListenerVector<RealType>& listener : listeners)
      listener.report(walker_index, name, ve_sample);
    for (const ListenerVector<RealType>& ion_listener : ion_listeners)
      ion_listener.report(walker_index, name, vi_sample);
    return value;
  };

  for (int iw = 0; iw < o_list.size(); iw++)
  {
    auto& coulomb_ab  = o_list.getCastedElement<CoulombPBCAB>(iw);
    coulomb_ab.value_ = evaluate_walker(iw, coulomb_ab, p_list[iw], listeners, ion_listeners);
  }
}

/** Evaluate the background term. Other constants are handled by AA potentials.
 *
 * \f$V_{bg}^{AB}=-\sum_{\alpha}\sum_{\beta} N^{\alpha} N^{\beta} q^{\alpha} q^{\beta} v_s(k=0) \f$
 * @todo Here is where the charge system has to be handled.
 */
CoulombPBCAB::Return_t CoulombPBCAB::evalConsts(const ParticleSet& Peln, bool report)
{
  int nelns = Peln.getTotalNum();
  int nions = pset_ions_.getTotalNum();
#if !defined(REMOVE_TRACEMANAGER)
  Ve_const.resize(nelns);
  Vi_const.resize(nions);
  Ve_const = 0.0;
  Vi_const = 0.0;
#endif
  mRealType Consts = 0.0;
  mRealType vs_k0  = AB->evaluateSR_k0();
  mRealType v1; //single particle energy
  for (int i = 0; i < nelns; ++i)
  {
    v1 = 0.0;
    for (int s = 0; s < NumSpeciesA; s++)
      v1 += NofSpeciesA[s] * Zspec[s];
    v1 *= -.5 * Qat[i] * vs_k0;
#if !defined(REMOVE_TRACEMANAGER)
    Ve_const(i) = v1;
#endif
    Consts += v1;
  }
  for (int i = 0; i < nions; ++i)
  {
    v1 = 0.0;
    for (int s = 0; s < NumSpeciesB; s++)
      v1 += NofSpeciesB[s] * Qspec[s];
    v1 *= -.5 * Zat[i] * vs_k0;
#if !defined(REMOVE_TRACEMANAGER)
    Vi_const(i) = v1;
#endif
    Consts += v1;
  }
  if (report)
    app_log() << "   Constant of PBCAB " << Consts << std::endl;
  return Consts;
}


CoulombPBCAB::Return_t CoulombPBCAB::evalSR(ParticleSet& P)
{
  constexpr mRealType czero(0);
  const auto& d_ab(P.getDistTableAB(myTableIndex));
  mRealType res = czero;
  //can be optimized but not important enough
  for (size_t b = 0; b < NptclB; ++b)
  {
    const auto& dist = d_ab.getDistRow(b);
    mRealType esum   = czero;
    for (size_t a = 0; a < NptclA; ++a)
      esum += Zat[a] * Vat[a]->splint(dist[a]) / dist[a];
    res += esum * Qat[b];
  }
  return res;
}


CoulombPBCAB::Return_t CoulombPBCAB::evalLR(ParticleSet& P)
{
  mRealType res = 0.0;
  const StructFact& RhoKA(pset_ions_.getSK());
  const StructFact& RhoKB(P.getSK());
  if (RhoKA.SuperCellEnum == SUPERCELL_SLAB)
    throw std::runtime_error(
        "CoulombPBCAB::evalLR RhoKA.SuperCellEnum == SUPERCELL_SLAB case not implemented. There was an implementation "
        "with complex-valued storage that may be resurrected using real-valued storage.");
  else
  {
    for (int i = 0; i < NumSpeciesA; i++)
    {
      mRealType esum = 0.0;
      for (int j = 0; j < NumSpeciesB; j++)
        esum += Qspec[j] *
            AB->evaluate(pset_ions_.getSimulationCell().getKLists().kshell, RhoKA.rhok_r[i], RhoKA.rhok_i[i],
                         RhoKB.rhok_r[j], RhoKB.rhok_i[j]);
      res += Zspec[i] * esum;
    }
  } //specion
  return res;
}


void CoulombPBCAB::initBreakup(ParticleSet& P)
{
  // This is ridiculous so many uneeded convenience variables that also make the object invalid at construction.
  // Illustraction of the terrible API and convoluted use of species attributes.
  SpeciesSet& tspeciesA(pset_ions_.getSpeciesSet());
  SpeciesSet& tspeciesB(P.getSpeciesSet());
  // Actually finding the index, code further on makes no sense if there wasn't already a charge attribute.
  int ChargeAttribIndxA = tspeciesA.addAttribute("charge");
  int ChargeAttribIndxB = tspeciesB.addAttribute("charge");
  NptclA                = pset_ions_.getTotalNum();
  NptclB                = P.getTotalNum();
  NumSpeciesA           = tspeciesA.TotalNum;
  NumSpeciesB           = tspeciesB.TotalNum;
  //Store information about charges and number of each species
  Zat.resize(NptclA);
  Zspec.resize(NumSpeciesA);
  Qat.resize(NptclB);
  Qspec.resize(NumSpeciesB);
  NofSpeciesA.resize(NumSpeciesA);
  NofSpeciesB.resize(NumSpeciesB);
  for (int spec = 0; spec < NumSpeciesA; spec++)
  {
    Zspec[spec]       = tspeciesA(ChargeAttribIndxA, spec);
    NofSpeciesA[spec] = pset_ions_.groupsize(spec);
  }
  for (int spec = 0; spec < NumSpeciesB; spec++)
  {
    Qspec[spec]       = tspeciesB(ChargeAttribIndxB, spec);
    NofSpeciesB[spec] = P.groupsize(spec);
  }
  for (int iat = 0; iat < NptclA; iat++)
    Zat[iat] = Zspec[pset_ions_.GroupID[iat]];
  for (int iat = 0; iat < NptclB; iat++)
    Qat[iat] = Qspec[P.GroupID[iat]];
  //    if(totQ>std::numeric_limits<RealType>::epsilon())
  //    {
  //      LOGMSG("PBCs not yet finished for non-neutral cells");
  //      OHMMS::Controller->abort();
  //    }
  ////Test if the box sizes are same (=> kcut same for fixed dimcut)
  kcdifferent =
      (std::abs(pset_ions_.getLattice().LR_kc - P.getLattice().LR_kc) > std::numeric_limits<RealType>::epsilon());
  minkc = std::min(pset_ions_.getLattice().LR_kc, P.getLattice().LR_kc);
  //AB->initBreakup(*PtclB);
  //initBreakup is called only once
  //AB = LRCoulombSingleton::getHandler(*PtclB);
  AB      = LRCoulombSingleton::getHandler(P);
  myConst = evalConsts(P);
  myRcut  = AB->get_rc(); //Basis.get_rc();
  // create the spline function for the short-range part assuming pure potential
  if (V0 == nullptr)
  {
    V0 = LRCoulombSingleton::createSpline4RbyVs(AB.get(), myRcut, myGrid.get());
    if (Vat.size())
    {
      APP_ABORT("CoulombPBCAB::initBreakup.  Vat is not empty\n");
    }
    Vat.resize(NptclA, V0.get());
    Vspec.resize(NumSpeciesA); //prepare for PP to overwrite it
  }

  //If ComputeForces is true, then we allocate space for the radial derivative functors.
  if (ComputeForces)
  {
    dAB = LRCoulombSingleton::getDerivHandler(P);
    if (fV0 == nullptr)
      fV0 = LRCoulombSingleton::createSpline4RbyVs(dAB.get(), myRcut, myGrid.get());
    if (dfV0 == nullptr)
      dfV0 = LRCoulombSingleton::createSpline4RbyVsDeriv(dAB.get(), myRcut, myGrid.get());
    if (fVat.size())
    {
      APP_ABORT("CoulombPBCAB::initBreakup.  Vat is not empty\n");
    }
    fVat.resize(NptclA, fV0.get());
    fdVat.resize(NptclA, dfV0.get());
    fVspec.resize(NumSpeciesA);
    fdVspec.resize(NumSpeciesA);
  }
}

/** add a local pseudo potential
 * @param groupID species index
 * @param ppot radial functor for \f$rV_{loc}\f$ on a grid
 */
void CoulombPBCAB::add(int groupID, std::unique_ptr<RadFunctorType>&& ppot)
{
  if (myGrid == nullptr)
  {
    myGrid = std::make_shared<LinearGrid<RealType>>();
    int ng = std::min(MaxGridPoints, static_cast<int>(myRcut / 1e-3) + 1);
    app_log() << "    CoulombPBCAB::add \n Setting a linear grid=[0," << myRcut << ") number of grid =" << ng
              << std::endl;
    myGrid->set(0, myRcut, ng);
  }
  if (Vspec[groupID] == nullptr)
  {
    V0.reset();

    app_log() << "    Creating the short-range pseudopotential for species " << groupID << std::endl;
    int ng = myGrid->size();
    std::vector<RealType> v(ng);
    for (int ig = 1; ig < ng - 2; ig++)
    {
      RealType r = (*myGrid)[ig];
      //need to multiply r for the LR
      v[ig] = -r * AB->evaluateLR(r) + ppot->splint(r);
    }
    v[0] = 2.0 * v[1] - v[2];
    //by construction, v has to go to zero at the boundary
    v[ng - 2]      = 0.0;
    v[ng - 1]      = 0.0;
    RealType deriv = (v[1] - v[0]) / ((*myGrid)[1] - (*myGrid)[0]);
    auto rfunc     = std::make_unique<RadFunctorType>(myGrid->makeClone(), v);
    rfunc->spline(0, deriv, ng - 1, 0.0);
    for (int iat = 0; iat < NptclA; iat++)
    {
      if (pset_ions_.GroupID[iat] == groupID)
        Vat[iat] = rfunc.get();
    }
    Vspec[groupID] = std::move(rfunc);
  }

  if (ComputeForces && fVspec[groupID] == nullptr)
  {
    app_log() << "    Creating the short-range pseudopotential derivatives for species " << groupID << std::endl;
    int ng = myGrid->size();
    //This is the value coming from optimized breakup for FORCES, not for energy.
    //Also.  Goal of this section is to create and store d/dr(rV), not d/dr(V)!!!
    std::vector<RealType> v(ng);
    std::vector<RealType> dv(ng);

    RealType ppot_val(0), ppot_deriv(0), ppot_2deriv(0);
    RealType lr_val(0), lr_deriv(0);
    for (int ig = 1; ig < ng - 2; ig++)
    {
      RealType r = (*myGrid)[ig];
      ppot_val   = ppot->splint(r, ppot_deriv, ppot_2deriv);
      lr_val     = dAB->evaluateLR(r);
      lr_deriv   = dAB->lrDf(r);

      v[ig]  = ppot_val - r * lr_val;
      dv[ig] = ppot_deriv - (lr_val + lr_deriv * r);
    }
    //This is a linear interpolation from v[2] and v[1] to v[0], assuming linear behavior.
    v[0]      = 2.0 * v[1] - v[2];
    v[ng - 2] = 0.0;
    v[ng - 1] = 0.0;

    dv[0]      = 2.0 * dv[1] - dv[2];
    dv[ng - 2] = 0;
    dv[ng - 1] = 0;

    auto ffunc  = std::make_unique<RadFunctorType>(myGrid->makeClone(), v);
    auto fdfunc = std::make_unique<RadFunctorType>(myGrid->makeClone(), dv);

    RealType fderiv = (dv[1] - dv[0]) / ((*myGrid)[1] - (*myGrid)[0]);

    ffunc->spline(0, dv[0], ng - 1, 0.0);
    fdfunc->spline(0, fderiv, ng - 1, 0.0);

    for (int iat = 0; iat < NptclA; iat++)
    {
      if (pset_ions_.GroupID[iat] == groupID)
      {
        fVat[iat]  = ffunc.get();
        fdVat[iat] = fdfunc.get();
      }
    }
    fVspec[groupID]  = std::move(ffunc);
    fdVspec[groupID] = std::move(fdfunc);
    //Done
  }

  if (ComputeForces)
  {
    if (OHMMS::Controller->rank() == 0)
    {
      FILE* fout = fopen("Vlocal.dat", "w");
      for (RealType r = 1.0e-8; r < myRcut; r += 1.0e-2)
      {
        RealType d_rV_dr, d2_rV_dr2;
        RealType Vr = Vat[0]->splint(r, d_rV_dr, d2_rV_dr2);
        Vr          = Vat[0]->splint(r);
        fprintf(fout, "%1.8e %1.12e %1.12e %1.12e\n", r, Vr, d_rV_dr, d2_rV_dr2);
      }
      fclose(fout);
    }
  }
}

CoulombPBCAB::Return_t CoulombPBCAB::evalLRwithForces(ParticleSet& P)
{
  //  const StructFact& RhoKA(pset_ions_.getSK());
  //  const StructFact& RhoKB(P.getSK());
  std::vector<TinyVector<RealType, DIM>> grad(pset_ions_.getTotalNum());
  for (int j = 0; j < NumSpeciesB; j++)
  {
    for (int iat = 0; iat < grad.size(); iat++)
      grad[iat] = TinyVector<RealType, DIM>(0.0, 0.0, 0.0);
    dAB->evaluateGrad(pset_ions_, P, j, Zat, grad);
    for (int iat = 0; iat < grad.size(); iat++)
      forces_[iat] += Qspec[j] * grad[iat];
  } // electron species
  return evalLR(P);
}

CoulombPBCAB::Return_t CoulombPBCAB::evalSRwithForces(ParticleSet& P)
{
  constexpr mRealType czero(0);
  const auto& d_ab(P.getDistTableAB(myTableIndex));
  mRealType res = czero;
  //Temporary variables for computing energy and forces.
  mRealType rV(0);
  mRealType frV(0), fdrV(0);
  mRealType rinv(1.0);
  //Magnitude of force.
  mRealType dvdr(0.0);
  for (size_t b = 0; b < NptclB; ++b)
  {
    const auto& dist = d_ab.getDistRow(b);
    const auto& dr   = d_ab.getDisplRow(b);
    mRealType esum   = czero;
    for (size_t a = 0; a < NptclA; ++a)
    {
      //Low hanging SIMD fruit here.  See J1/J2 grad computation.
      rinv = 1.0 / dist[a];
      rV   = Vat[a]->splint(dist[a]);
      frV  = fVat[a]->splint(dist[a]);
      fdrV = fdVat[a]->splint(dist[a]);
      dvdr = Qat[b] * Zat[a] * (fdrV - frV * rinv) * rinv;
      forces_[a][0] -= dvdr * dr[a][0] * rinv;
      forces_[a][1] -= dvdr * dr[a][1] * rinv;
      forces_[a][2] -= dvdr * dr[a][2] * rinv;
      esum += Zat[a] * rV * rinv; //Total energy accumulation
    }
    res += esum * Qat[b];
  }
  return res;
}

void CoulombPBCAB::evalPerParticleConsts(Vector<RealType>& pp_consts_src, Vector<RealType>& pp_consts_trg) const
{
  int nelns = Peln.getTotalNum();
  int nions = pset_ions_.getTotalNum();
  pp_consts_trg.resize(nelns, 0.0);
  pp_consts_src.resize(nions, 0.0);

  RealType vs_k0 = AB->evaluateSR_k0();
  RealType v1; //single particle energy
  for (int i = 0; i < nelns; ++i)
  {
    v1 = 0.0;
    for (int s = 0; s < NumSpeciesA; s++)
      v1 += NofSpeciesA[s] * Zspec[s];
    v1 *= -.5 * Qat[i] * vs_k0;
    pp_consts_trg[i] = v1;
  }
  for (int i = 0; i < nions; ++i)
  {
    v1 = 0.0;
    for (int s = 0; s < NumSpeciesB; s++)
      v1 += NofSpeciesB[s] * Qspec[s];
    v1 *= -.5 * Zat[i] * vs_k0;
    pp_consts_src[i] = v1;
  }
}

void CoulombPBCAB::createResource(ResourceCollection& collection) const
{
  auto new_res = std::make_unique<CoulombPBCABMultiWalkerResource>();
  evalPerParticleConsts(new_res->pp_consts_src, new_res->pp_consts_trg);
  auto resource_index = collection.addResource(std::move(new_res));
}

void CoulombPBCAB::acquireResource(ResourceCollection& collection,
                                   const RefVectorWithLeader<OperatorBase>& o_list) const
{
  auto& o_leader          = o_list.getCastedLeader<CoulombPBCAB>();
  o_leader.mw_res_handle_ = collection.lendResource<CoulombPBCABMultiWalkerResource>();
}

void CoulombPBCAB::releaseResource(ResourceCollection& collection,
                                   const RefVectorWithLeader<OperatorBase>& o_list) const
{
  auto& o_leader = o_list.getCastedLeader<CoulombPBCAB>();
  collection.takebackResource(o_leader.mw_res_handle_);
}

} // namespace qmcplusplus
