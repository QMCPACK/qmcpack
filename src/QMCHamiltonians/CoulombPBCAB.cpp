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


#include "CoulombPBCAB.h"
#include "Particle/DistanceTableData.h"
#include "Message/Communicate.h"
#include "Utilities/ProgressReportEngine.h"

namespace qmcplusplus
{
CoulombPBCAB::CoulombPBCAB(ParticleSet& ions, ParticleSet& elns, bool computeForces)
    : ForceBase(ions, elns),
      PtclA(ions),
      myTableIndex(elns.addTable(ions)),
      myConst(0.0),
      myGrid(nullptr),
      V0(nullptr),
      fV0(nullptr),
      dfV0(nullptr),
      ComputeForces(computeForces),
      MaxGridPoints(10000),
      Peln(elns),
      Pion(ions)
{
  ReportEngine PRE("CoulombPBCAB", "CoulombPBCAB");
  set_energy_domain(potential);
  two_body_quantum_domain(ions, elns);
  if (ComputeForces)
    PtclA.turnOnPerParticleSK();
  initBreakup(elns);
  prefix = "Flocal";
  app_log() << "  Rcut                " << myRcut << std::endl;
  app_log() << "  Maximum K shell     " << AB->MaxKshell << std::endl;
  app_log() << "  Number of k vectors " << AB->Fk.size() << std::endl;
}

OperatorBase* CoulombPBCAB::makeClone(ParticleSet& qp, TrialWaveFunction& psi)
{
  CoulombPBCAB* myclone    = new CoulombPBCAB(PtclA, qp, ComputeForces);
  myclone->FirstForceIndex = FirstForceIndex;
  if (myGrid)
    myclone->myGrid = new GridType(*myGrid);
  for (int ig = 0; ig < Vspec.size(); ++ig)
  {
    if (Vspec[ig])
    {
      RadFunctorType* apot = Vspec[ig]->makeClone();
      myclone->Vspec[ig]   = apot;
      for (int iat = 0; iat < PtclA.getTotalNum(); ++iat)
      {
        if (PtclA.GroupID[iat] == ig)
          myclone->Vat[iat] = apot;
      }
    }
  }
  //If forces exist, force arrays will have been allocated.  Iterate over one
  //such array to clone.
  for (int ig = 0; ig < fVspec.size(); ig++)
  {
    if (fVspec[ig])
    {
      RadFunctorType* apot  = fVspec[ig]->makeClone();
      RadFunctorType* dapot = fdVspec[ig]->makeClone();
      myclone->fVspec[ig]   = apot;
      myclone->fdVspec[ig]  = dapot;
      for (int iat = 0; iat < PtclA.getTotalNum(); ++iat)
      {
        if (PtclA.GroupID[iat] == ig)
        {
          myclone->fVat[iat]  = apot;
          myclone->fdVat[iat] = dapot;
        }
      }
    }
  }
  return myclone;
}

CoulombPBCAB::~CoulombPBCAB()
{
  delete V0;
  delete fV0;
  delete dfV0;

  V0   = nullptr;
  fV0  = nullptr;
  dfV0 = nullptr;

  //This takes care of species dependent terms.
  for (int ig = 0; ig < Vspec.size(); ig++)
  {
    delete Vspec[ig];
    Vspec[ig] = nullptr;
  }
  //If forces were initialized, fVspec.size()!=0.  This cleans up.
  for (int ig = 0; ig < fVspec.size(); ig++)
  {
    delete fVspec[ig];
    delete fdVspec[ig];
    fVspec[ig]  = nullptr;
    fdVspec[ig] = nullptr;
  }
}

void CoulombPBCAB::resetTargetParticleSet(ParticleSet& P)
{
  int tid = P.addTable(PtclA);
  if (tid != myTableIndex)
  {
    APP_ABORT("CoulombPBCAB::resetTargetParticleSet found inconsistent table index");
  }
  AB->resetTargetParticleSet(P);
}

void CoulombPBCAB::addObservables(PropertySetType& plist, BufferType& collectables)
{
  myIndex = plist.add(myName.c_str());
  if (ComputeForces)
    addObservablesF(plist);
}


#if !defined(REMOVE_TRACEMANAGER)
void CoulombPBCAB::contribute_particle_quantities() { request.contribute_array(myName); }

void CoulombPBCAB::checkout_particle_quantities(TraceManager& tm)
{
  streaming_particles = request.streaming_array(myName);
  if (streaming_particles)
  {
    Ve_sample = tm.checkout_real<1>(myName, Peln);
    Vi_sample = tm.checkout_real<1>(myName, Pion);
  }
}

void CoulombPBCAB::delete_particle_quantities()
{
  if (streaming_particles)
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
    forces = 0.0;
    Value  = evalLRwithForces(P) + evalSRwithForces(P) + myConst;
  }
  else
#if !defined(REMOVE_TRACEMANAGER)
      if (streaming_particles)
    Value = evaluate_sp(P);
  else
#endif
    Value = evalLR(P) + evalSR(P) + myConst;
  return Value;
}

CoulombPBCAB::Return_t CoulombPBCAB::evaluateWithIonDerivs(ParticleSet& P,
                                                           ParticleSet& ions,
                                                           TrialWaveFunction& psi,
                                                           ParticleSet::ParticlePos_t& hf_terms,
                                                           ParticleSet::ParticlePos_t& pulay_terms)
{
  if (ComputeForces)
  {
    forces = 0.0;
    Value  = evalLRwithForces(P) + evalSRwithForces(P) + myConst;
    hf_terms -= forces;
    //And no Pulay contribution.
  }
  else
    Value = evalLR(P) + evalSR(P) + myConst;
  return Value;
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
    const DistanceTableData& d_ab(P.getDistTable(myTableIndex));
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
    const StructFact& RhoKA(*(PtclA.SK));
    const StructFact& RhoKB(*(P.SK));
    if (RhoKA.SuperCellEnum == SUPERCELL_SLAB)
    {
      APP_ABORT("CoulombPBCAB::evaluate_sp single particle traces have not been implemented for slab geometry");
    }
    else
    {
      //jtk mark: needs optimizations for USE_REAL_STRUCT_FACTOR
      //          will likely require new function definitions
      RealType v1; //single particle energy
      RealType q;
      for (int i = 0; i < P.getTotalNum(); ++i)
      {
        q  = .5 * Qat[i];
        v1 = 0.0;
        for (int s = 0; s < NumSpeciesA; s++)
#if defined(USE_REAL_STRUCT_FACTOR)
          v1 += Zspec[s] * q *
              AB->evaluate(RhoKA.KLists.kshell, RhoKA.rhok_r[s], RhoKA.rhok_i[s], RhoKB.eikr_r[i], RhoKB.eikr_i[i]);
#else
          v1 += Zspec[s] * q * AB->evaluate(RhoKA.KLists.kshell, RhoKA.rhok[s], RhoKB.eikr[i]);
#endif
        Ve_samp(i) += v1;
        Vlr += v1;
      }
      for (int i = 0; i < PtclA.getTotalNum(); ++i)
      {
        q  = .5 * Zat[i];
        v1 = 0.0;
        for (int s = 0; s < NumSpeciesB; s++)
#if defined(USE_REAL_STRUCT_FACTOR)
          v1 += Qspec[s] * q *
              AB->evaluate(RhoKB.KLists.kshell, RhoKB.rhok_r[s], RhoKB.rhok_i[s], RhoKA.eikr_r[i], RhoKA.eikr_i[i]);
#else
          v1 += Qspec[s] * q * AB->evaluate(RhoKB.KLists.kshell, RhoKB.rhok[s], RhoKA.eikr[i]);
#endif
        Vi_samp(i) += v1;
        Vlr += v1;
      }
    }
  }
  for (int i = 0; i < Ve_samp.size(); ++i)
    Ve_samp(i) += Ve_const(i);
  for (int i = 0; i < Vi_samp.size(); ++i)
    Vi_samp(i) += Vi_const(i);
  Value = Vsr + Vlr + Vc;
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
  return Value;
}
#endif


/** Evaluate the background term. Other constants are handled by AA potentials.
 *
 * \f$V_{bg}^{AB}=-\sum_{\alpha}\sum_{\beta} N^{\alpha} N^{\beta} q^{\alpha} q^{\beta} v_s(k=0) \f$
 * @todo Here is where the charge system has to be handled.
 */
CoulombPBCAB::Return_t CoulombPBCAB::evalConsts(bool report)
{
  int nelns = Peln.getTotalNum();
  int nions = Pion.getTotalNum();
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
  const DistanceTableData& d_ab(P.getDistTable(myTableIndex));
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
  const StructFact& RhoKA(*(PtclA.SK));
  const StructFact& RhoKB(*(P.SK));
  if (RhoKA.SuperCellEnum == SUPERCELL_SLAB)
  {
    const DistanceTableData& d_ab(P.getDistTable(myTableIndex));
    for (int iat = 0; iat < NptclA; ++iat)
    {
      mRealType u = 0;
#if !defined(USE_REAL_STRUCT_FACTOR)
      const int slab_dir = OHMMS_DIM - 1;
      for (int nn = d_ab.M[iat], jat = 0; nn < d_ab.M[iat + 1]; ++nn, ++jat)
        u += Qat[jat] * AB->evaluate_slab(d_ab.dr(nn)[slab_dir], RhoKA.KLists.kshell, RhoKA.eikr[iat], RhoKB.eikr[jat]);
#endif
      res += Zat[iat] * u;
    }
  }
  else
  {
    for (int i = 0; i < NumSpeciesA; i++)
    {
      mRealType esum = 0.0;
      for (int j = 0; j < NumSpeciesB; j++)
      {
#if defined(USE_REAL_STRUCT_FACTOR)
        esum += Qspec[j] *
            AB->evaluate(RhoKA.KLists.kshell, RhoKA.rhok_r[i], RhoKA.rhok_i[i], RhoKB.rhok_r[j], RhoKB.rhok_i[j]);
#else
        esum += Qspec[j] * AB->evaluate(RhoKA.KLists.kshell, RhoKA.rhok[i], RhoKB.rhok[j]);
#endif
      } //speceln
      res += Zspec[i] * esum;
    }
  } //specion
  return res;
}


void CoulombPBCAB::initBreakup(ParticleSet& P)
{
  SpeciesSet& tspeciesA(PtclA.getSpeciesSet());
  SpeciesSet& tspeciesB(P.getSpeciesSet());
  int ChargeAttribIndxA = tspeciesA.addAttribute("charge");
  int MemberAttribIndxA = tspeciesA.addAttribute("membersize");
  int ChargeAttribIndxB = tspeciesB.addAttribute("charge");
  int MemberAttribIndxB = tspeciesB.addAttribute("membersize");
  NptclA                = PtclA.getTotalNum();
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
    NofSpeciesA[spec] = static_cast<int>(tspeciesA(MemberAttribIndxA, spec));
  }
  for (int spec = 0; spec < NumSpeciesB; spec++)
  {
    Qspec[spec]       = tspeciesB(ChargeAttribIndxB, spec);
    NofSpeciesB[spec] = static_cast<int>(tspeciesB(MemberAttribIndxB, spec));
  }
  RealType totQ = 0.0;
  for (int iat = 0; iat < NptclA; iat++)
    totQ += Zat[iat] = Zspec[PtclA.GroupID[iat]];
  for (int iat = 0; iat < NptclB; iat++)
    totQ += Qat[iat] = Qspec[P.GroupID[iat]];
  //    if(totQ>std::numeric_limits<RealType>::epsilon())
  //    {
  //      LOGMSG("PBCs not yet finished for non-neutral cells");
  //      OHMMS::Controller->abort();
  //    }
  ////Test if the box sizes are same (=> kcut same for fixed dimcut)
  kcdifferent = (std::abs(PtclA.Lattice.LR_kc - P.Lattice.LR_kc) > std::numeric_limits<RealType>::epsilon());
  minkc       = std::min(PtclA.Lattice.LR_kc, P.Lattice.LR_kc);
  //AB->initBreakup(*PtclB);
  //initBreakup is called only once
  //AB = LRCoulombSingleton::getHandler(*PtclB);
  AB      = LRCoulombSingleton::getHandler(P);
  myConst = evalConsts();
  myRcut  = AB->get_rc(); //Basis.get_rc();
  // create the spline function for the short-range part assuming pure potential
  if (V0 == nullptr)
  {
    V0 = LRCoulombSingleton::createSpline4RbyVs(AB.get(), myRcut, myGrid);
    if (Vat.size())
    {
      APP_ABORT("CoulombPBCAB::initBreakup.  Vat is not empty\n");
    }
    Vat.resize(NptclA, V0);
    Vspec.resize(NumSpeciesA, nullptr); //prepare for PP to overwrite it
  }

  //If ComputeForces is true, then we allocate space for the radial derivative functors.
  if (ComputeForces)
  {
    dAB = LRCoulombSingleton::getDerivHandler(P);
    if (fV0 == nullptr)
      fV0 = LRCoulombSingleton::createSpline4RbyVs(dAB.get(), myRcut, myGrid);
    if (dfV0 == nullptr)
      dfV0 = LRCoulombSingleton::createSpline4RbyVsDeriv(dAB.get(), myRcut, myGrid);
    if (fVat.size())
    {
      APP_ABORT("CoulombPBCAB::initBreakup.  Vat is not empty\n");
    }
    fVat.resize(NptclA, fV0);
    fdVat.resize(NptclA, dfV0);
    fVspec.resize(NumSpeciesA, nullptr);
    fdVspec.resize(NumSpeciesA, nullptr);
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
    myGrid = new LinearGrid<RealType>;
    int ng = std::min(MaxGridPoints, static_cast<int>(myRcut / 1e-3) + 1);
    app_log() << "    CoulombPBCAB::add \n Setting a linear grid=[0," << myRcut << ") number of grid =" << ng
              << std::endl;
    myGrid->set(0, myRcut, ng);
  }
  if (Vspec[groupID] == nullptr)
  {
    delete V0;
    V0 = nullptr;

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
    v[ng - 2]             = 0.0;
    v[ng - 1]             = 0.0;
    RadFunctorType* rfunc = new RadFunctorType(myGrid, v);
    RealType deriv        = (v[1] - v[0]) / ((*myGrid)[1] - (*myGrid)[0]);
    rfunc->spline(0, deriv, ng - 1, 0.0);
    Vspec[groupID] = rfunc;
    for (int iat = 0; iat < NptclA; iat++)
    {
      if (PtclA.GroupID[iat] == groupID)
        Vat[iat] = rfunc;
    }
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

    RadFunctorType* ffunc  = new RadFunctorType(myGrid, v);
    RadFunctorType* fdfunc = new RadFunctorType(myGrid, dv);

    RealType fderiv = (dv[1] - dv[0]) / ((*myGrid)[1] - (*myGrid)[0]);

    ffunc->spline(0, dv[0], ng - 1, 0.0);
    fdfunc->spline(0, fderiv, ng - 1, 0.0);

    fVspec[groupID]  = ffunc;
    fdVspec[groupID] = fdfunc;
    for (int iat = 0; iat < NptclA; iat++)
    {
      if (PtclA.GroupID[iat] == groupID)
      {
        fVat[iat]  = ffunc;
        fdVat[iat] = fdfunc;
      }
    }
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
  //  const StructFact& RhoKA(*(PtclA.SK));
  //  const StructFact& RhoKB(*(P.SK));
  std::vector<TinyVector<RealType, DIM>> grad(PtclA.getTotalNum());
  for (int j = 0; j < NumSpeciesB; j++)
  {
    for (int iat = 0; iat < grad.size(); iat++)
      grad[iat] = TinyVector<RealType, DIM>(0.0, 0.0, 0.0);
    dAB->evaluateGrad(PtclA, P, j, Zat, grad);
    for (int iat = 0; iat < grad.size(); iat++)
      forces[iat] += Qspec[j] * grad[iat];
  } // electron species
  return evalLR(P);
}

CoulombPBCAB::Return_t CoulombPBCAB::evalSRwithForces(ParticleSet& P)
{
  constexpr mRealType czero(0);
  const DistanceTableData& d_ab(P.getDistTable(myTableIndex));
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
      forces[a][0] -= dvdr * dr[a][0] * rinv;
      forces[a][1] -= dvdr * dr[a][1] * rinv;
      forces[a][2] -= dvdr * dr[a][2] * rinv;
      esum += Zat[a] * rV * rinv; //Total energy accumulation
    }
    res += esum * Qat[b];
  }
  return res;
}
} // namespace qmcplusplus
