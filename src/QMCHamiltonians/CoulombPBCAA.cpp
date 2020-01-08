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


#include "QMCHamiltonians/CoulombPBCAA.h"
#include "Particle/DistanceTableData.h"
#include "Utilities/ProgressReportEngine.h"
#include <numeric>

namespace qmcplusplus
{


CoulombPBCAA::CoulombPBCAA(ParticleSet& ref, bool active, bool computeForces)
    : ForceBase(ref, ref),
      AA(0),
      myGrid(0),
      rVs(0),
      dAA(0),
      myGridforce(0),
      rVsforce(0),
      is_active(active),
      FirstTime(true),
      myConst(0.0),
      ComputeForces(computeForces),
      Ps(ref),
      d_aa_ID(ref.addTable(ref, DT_SOA_PREFERRED))
{
  ReportEngine PRE("CoulombPBCAA", "CoulombPBCAA");
  set_energy_domain(potential);
  two_body_quantum_domain(ref);
  PtclRefName = ref.getDistTable(d_aa_ID).getName();
  initBreakup(ref);

  if (ComputeForces)
  {
    ref.turnOnPerParticleSK();
    update_source(ref);
  }


  // Setup anistropic ewald
  ewaldref::RealMat A;
  ewaldref::ChargeArray Q;

  A = Ps.Lattice.R;

  Q.resize(NumCenters);
  for(int i=0;i<NumCenters;++i)
    Q[i] = Zat[i];



  //RealType qmcpack_kappa = AA->LR_kc;
  RealType qmcpack_sigma = std::sqrt(AA->LR_kc / (2.0 * AA->LR_rc));
  RealType qmcpack_kappa = 1./(std::sqrt(2.0)*qmcpack_sigma);
  ewaldtools::AnisotropicEwald ewald_qp(A,Q,1e-10,qmcpack_kappa);

  //ewald.initialize(A,Q);
  ewald.initialize(A,Q,1e-10,qmcpack_kappa);


  ewaldref::IntVec nmax = 20;
  ewald.setNmax(nmax);

  if (!is_active)
  {
    update_source(ref);

    ewaldref::PosArray R;
    R.resize(NumCenters);
    for(int i=0;i<NumCenters;++i)
      R[i] = Ps.R[i];
    
    RealType eE   = ewald.ewaldEnergy(R);
    RealType eEf  = ewald.fixedGridEwaldEnergy(R);
    //RealType eEt  = ewald.ewaldEnergy(Ps.R);
    //RealType eEft = ewald.fixedGridEwaldEnergy(Ps.R);

    RealType c  = ewald.ewaldEnergyConst();
    RealType sr = ewald.ewaldEnergySR(Ps.R);
    RealType lr = ewald.ewaldEnergyLR(Ps.R);

    const DistanceTableData& dt(Ps.getDistTable(d_aa_ID));
    //RealType srdt = ewald.ewaldEnergySRDT(dt);
    //RealType lrdt = ewald.ewaldEnergyLRDT(dt);


    //ewald_qp.madelungEnergy0();
    RealType eE_qp = ewald_qp.ewaldEnergy(Ps.R);
    //RealType c_qp   = ewald_qp.ewaldEnergyConst();
    //RealType sr_qp  = ewald_qp.ewaldEnergySRDT(dt);
    //RealType lr_qp  = ewald_qp.ewaldEnergyLRDT(dt);
    //RealType sr0_qp = ewald_qp.ewaldEnergySR0DT(dt);

    RealType c_qmc = myConst;
    RealType sr_qmc = evalSR(Ps);
    RealType lr_qmc = evalLR(Ps);
    RealType E_qmc  = c_qmc + sr_qmc + lr_qmc;


    app_log()<<std::setprecision(14);
    app_log()<<std::endl;
    app_log()<<"   adaptive aniso ewald energy: "<<eE<<std::endl;
    app_log()<<"   fixed aniso ewald energy   : "<<eEf<<std::endl;

    //app_log()<<"   adaptive aniso ewald energy: "<<eEt<<std::endl;
    //app_log()<<"   fixed aniso ewald energy   : "<<eEft<<std::endl;
    app_log()<<"   ewald energy c+sr+lr       : "<<c+sr+lr<<std::endl;
    //app_log()<<"   ewald energy c+sr+lr dt    : "<<c+srdt+lrdt<<std::endl;
    app_log()<<"   aniso ewald nmax           : "<<ewald.getNmax()<<std::endl;
    app_log()<<std::setprecision(14);
    //app_log()<<"   qmcpack constant : "<<c_qmc<<std::endl;
    //app_log()<<"   ewald   constant : "<<c_qp<<std::endl;
    //app_log()<<"   qmcpack SR  : "<<sr_qmc<<std::endl;
    //app_log()<<"   ewald   SR  : "<<sr_qp<<std::endl;
    //app_log()<<"   ewald   SR0 : "<<sr0_qp<<std::endl;
    //app_log()<<"   qmcpack LR  : "<<lr_qmc<<std::endl;
    //app_log()<<"   ewald   LR  : "<<lr_qp<<std::endl;
    //app_log()<<"   qmcpack tot : "<<c_qmc+sr_qmc+lr_qmc<<std::endl;
    //app_log()<<"   ewald   tot : "<< c_qp+sr_qp+lr_qp <<std::endl;
    //app_log()<<"   ewald   tot0: "<< c_qp+sr0_qp+lr_qp <<std::endl;
    //app_log()<<std::endl;
    //app_log()<<"   SR error    : "<<sr_qmc-sr_qp<<std::endl;
    //app_log()<<"   tot error   : "<<c_qmc+sr_qmc+lr_qmc-c_qp-sr_qp-lr_qp<<std::endl;
    //app_log()<<std::endl;
    app_log()<<"   ewald_qp aniso energy      : "<<eE<<std::endl;
    app_log()<<"   aniso ewald nmax           : "<<ewald_qp.getNmax()<<std::endl;
    app_log()<<std::setprecision(14);

    
    //ewaldtools::IntVec nmax = 20;
    //ewald.setupOpt(nmax);

    //auto& ewald_opt = ewald_qp;
    auto& ewald_opt = ewald;

    RealType c_ref  = ewald_opt.ewaldEnergyConst();
    RealType sr_ref = ewald_opt.ewaldEnergySR(Ps.R);
    RealType lr_ref = ewald_opt.ewaldEnergyLR(Ps.R);
    //ewald_opt.setupOpt();
    ewald_opt.setupOpt(nmax);
    RealType c_opt  = ewald_opt.ewaldEnergyConst();
    RealType sr_opt = ewald_opt.ewaldEnergySROpt(dt);
    RealType lr_opt = ewald_opt.ewaldEnergyLROpt(Ps.R);
    app_log()<<"   ewald opt energy           : "<<c_opt+sr_opt+lr_opt<<std::endl;
    app_log()<<"   ewald opt energy2          : "<<ewald.ewaldEnergyOpt(Ps.R,dt)<<std::endl;
    app_log()<<"   qmcpack   energy           : "<<E_qmc<<std::endl;
    app_log()<<"   ewald     LR  : "<<lr_ref<<std::endl;
    app_log()<<"   ewald opt LR  : "<<lr_opt<<std::endl;
    app_log()<<"   ewald     SR  : "<<sr_ref<<std::endl;
    app_log()<<"   ewald opt SR  : "<<sr_opt<<std::endl;
    app_log()<<"   qmcpack   SR  : "<<sr_qmc<<std::endl;

    app_log()<<std::endl;
    app_log()<<std::endl;

    RealType Vii_ref        = ewaldref::ewaldEnergy(A, R, Q);
    RealType Vdiff_per_atom = std::abs(Value - Vii_ref) / NumCenters;
    app_log() << "Checking ion-ion Ewald energy against reference..." << std::endl;
    if (Vdiff_per_atom > Ps.Lattice.LR_tol)
    {
      app_log() << std::setprecision(14);
      app_log() << std::endl;
      app_log() << "Error in ion-ion Ewald energy exceeds " << Ps.Lattice.LR_tol << " Ha/atom tolerance." << std::endl;
      app_log() << std::endl;
      app_log() << "  Reference ion-ion energy: " << Vii_ref << std::endl;
      app_log() << "  QMCPACK   ion-ion energy: " << Value << std::endl;
      app_log() << "            ion-ion diff  : " << Value - Vii_ref << std::endl;
      app_log() << "            diff/atom     : " << (Value - Vii_ref) / NumCenters << std::endl;
      app_log() << "            tolerance     : " << Ps.Lattice.LR_tol << std::endl;
      app_log() << std::endl;
      app_log() << "Please try increasing the LR_dim_cutoff parameter in the <simulationcell/>" << std::endl;
      app_log() << "input.  Alternatively, the tolerance can be increased by setting the" << std::endl;
      app_log() << "LR_tol parameter in <simulationcell/> to a value greater than " << Ps.Lattice.LR_tol << ". "
                << std::endl;
      app_log() << "If you increase the tolerance, please perform careful checks of energy" << std::endl;
      app_log() << "differences to ensure this error is controlled for your application." << std::endl;
      app_log() << std::endl;

      APP_ABORT("ion-ion check failed")
    }
    else
    {
      app_log() << "  Check passed." << std::endl;
    }

    //APP_ABORT("explore")

  }
  prefix = "F_AA";
  app_log() << "  Maximum K shell " << AA->MaxKshell << std::endl;
  app_log() << "  Number of k vectors " << AA->Fk.size() << std::endl;
  app_log() << "  Fixed Coulomb potential for " << ref.getName();
  app_log() << "\n    e-e Madelung Const. =" << MC0 << "\n    Vtot     =" << Value << std::endl;
}

CoulombPBCAA::~CoulombPBCAA() {}

void CoulombPBCAA::addObservables(PropertySetType& plist, BufferType& collectables)
{
  addValue(plist);
  if (ComputeForces)
    addObservablesF(plist);
}

void CoulombPBCAA::update_source(ParticleSet& s)
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
  NewValue = Value = eL + eS + myConst;
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
void CoulombPBCAA::contribute_particle_quantities() { request.contribute_array(myName); }

void CoulombPBCAA::checkout_particle_quantities(TraceManager& tm)
{
  streaming_particles = request.streaming_array(myName);
  if (streaming_particles)
  {
    V_sample = tm.checkout_real<1>(myName, Ps);
    if (!is_active)
      evaluate_sp(Ps);
  }
}

void CoulombPBCAA::delete_particle_quantities()
{
  if (streaming_particles)
    delete V_sample;
}
#endif


CoulombPBCAA::Return_t CoulombPBCAA::evaluate(ParticleSet& P)
{
  ScopedTimer evaltimer(TimerManager.createTimer("CoulombPBCAA_eval_total"));

  const DistanceTableData& dt(P.getDistTable(d_aa_ID));

  if (is_active)
  {
    {
    ScopedTimer evaltimer(TimerManager.createTimer("QMCPACK eval"));

#if !defined(REMOVE_TRACEMANAGER)
    if (streaming_particles)
      Value = evaluate_sp(P);
    else
#endif
      Value = evalLR(P) + evalSR(P) + myConst;
    }

    RealType eval;
    {
      ScopedTimer evaltimer(TimerManager.createTimer("Aniso eval"));
      //eval = ewald.fixedGridEwaldEnergy(P.R);
      eval = ewald.ewaldEnergyOpt(P.R,dt);
    }

    app_log()<<std::endl;
    app_log()<<"  QMCPACK value: "<<Value<<std::endl;
    app_log()<<"  Aniso   value: "<<eval<<std::endl;
    app_log()<<"  QMCPACK    SR: "<<evalSR(P)<<std::endl;
    app_log()<<"  Aniso      SR: "<<ewald.ewaldEnergySROpt(dt)<<std::endl;
    app_log()<<"  QMCPACK    LR: "<<evalLR(P)+myConst<<std::endl;
    app_log()<<"  Aniso      LR: "<<ewald.ewaldEnergyLROpt(P.R)<<std::endl;
    app_log()<<std::endl;
  }
  return Value;
}

CoulombPBCAA::Return_t CoulombPBCAA::evaluateWithIonDerivs(ParticleSet& P,
                                                           ParticleSet& ions,
                                                           TrialWaveFunction& psi,
                                                           ParticleSet::ParticlePos_t& hf_terms,
                                                           ParticleSet::ParticlePos_t& pulay_terms)
{
  hf_terms -= forces;
  //No pulay term.
  return Value;
}

#if !defined(REMOVE_TRACEMANAGER)
CoulombPBCAA::Return_t CoulombPBCAA::evaluate_sp(ParticleSet& P)
{
  mRealType Vsr              = 0.0;
  mRealType Vlr              = 0.0;
  mRealType& Vc              = myConst;
  Array<RealType, 1>& V_samp = V_samp_tmp;
  V_samp                     = 0.0;
  {
    //SR
    const DistanceTableData& d_aa(P.getDistTable(d_aa_ID));
    RealType z;
    if (d_aa.DTType == DT_SOA)
    {
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
    }
    else
    {
#ifndef ENABLE_SOA
      for (int ipart = 0; ipart < NumCenters; ipart++)
      {
        z = .5 * Zat[ipart];
        for (int nn = d_aa.M[ipart], jpart = ipart + 1; nn < d_aa.M[ipart + 1]; nn++, jpart++)
        {
          RealType pairpot = z * Zat[jpart] * d_aa.rinv(nn) * rVs->splint(d_aa.r(nn));
          V_samp(ipart) += pairpot;
          V_samp(jpart) += pairpot;
          Vsr += pairpot;
        }
      }
#endif
    }
    Vsr *= 2.0;
  }
  {
    //LR
    const StructFact& PtclRhoK(*(P.SK));
    if (PtclRhoK.SuperCellEnum == SUPERCELL_SLAB)
    {
      APP_ABORT("CoulombPBCAA::evaluate_sp single particle traces have not been implemented for slab geometry");
    }
    else
    {
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
              AA->evaluate(PtclRhoK.KLists.kshell, PtclRhoK.rhok_r[s], PtclRhoK.rhok_i[s], PtclRhoK.eikr_r[i],
                           PtclRhoK.eikr_i[i]);
#else
          v1 += z * Zspec[s] * AA->evaluate(PtclRhoK.KLists.kshell, PtclRhoK.rhok[s], PtclRhoK.eikr[i]);
#endif
        }
        V_samp(i) += v1;
        Vlr += v1;
      }
    }
  }
  for (int i = 0; i < V_samp.size(); ++i)
    V_samp(i) += V_const(i);
  Value = Vsr + Vlr + Vc;
#if defined(TRACE_CHECK)
  RealType Vlrnow = evalLR(P);
  RealType Vsrnow = evalSR(P);
  RealType Vcnow  = myConst;
  RealType Vnow   = Vlrnow + Vsrnow + Vcnow;
  RealType Vsum   = V_samp.sum();
  RealType Vcsum  = V_const.sum();
  RealType Vsrold = evalSR_old(P);
  RealType Vlrold = evalLR_old(P);
  RealType Vcold  = evalConsts_old(false);
  RealType Vcorig = evalConsts_orig(false);
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
  if (std::abs(Vsrold - Vsrnow) > TraceManager::trace_tol)
  {
    app_log() << "versiontest: CoulombPBCAA::evalSR()" << std::endl;
    app_log() << "versiontest:    old:" << Vsrold << std::endl;
    app_log() << "versiontest:    mod:" << Vsrnow << std::endl;
    APP_ABORT("Trace check failed");
  }
  if (std::abs(Vlrold - Vlrnow) > TraceManager::trace_tol)
  {
    app_log() << "versiontest: CoulombPBCAA::evalLR()" << std::endl;
    app_log() << "versiontest:    old:" << Vlrold << std::endl;
    app_log() << "versiontest:    mod:" << Vlrnow << std::endl;
    APP_ABORT("Trace check failed");
  }
  if (std::abs(Vcold - Vcorig) > TraceManager::trace_tol || std::abs(Vcnow - Vcorig) > TraceManager::trace_tol)
  {
    app_log() << "versiontest: CoulombPBCAA::evalConsts()" << std::endl;
    app_log() << "versiontest:    old:" << Vcold << std::endl;
    app_log() << "versiontest:   orig:" << Vcorig << std::endl;
    app_log() << "versiontest:    mod:" << Vcnow << std::endl;
    APP_ABORT("Trace check failed");
  }
#endif
  return Value;
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
  V_samp_tmp.resize(NumCenters);
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
  if (rVs == 0)
  {
    rVs = LRCoulombSingleton::createSpline4RbyVs(AA, myRcut, myGrid);
  }
  if (ComputeForces)
  {
    dAA = LRCoulombSingleton::getDerivHandler(P);
    if (rVsforce == 0)
    {
      rVsforce = LRCoulombSingleton::createSpline4RbyVs(dAA, myRcut, myGridforce);
    }
  }

  P.update();
}


CoulombPBCAA::Return_t CoulombPBCAA::evalLRwithForces(ParticleSet& P)
{
  //  const StructFact& PtclRhoK(*(P.SK));
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
  const DistanceTableData& d_aa(P.getDistTable(d_aa_ID));
  mRealType SR = 0.0;
  if (d_aa.DTType == DT_SOA)
  {
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
  }
  else
  {
#ifndef ENABLE_SOA
    for (int ipart = 0; ipart < NumCenters; ipart++)
    {
      mRealType esum = 0.0;
      for (int nn = d_aa.M[ipart], jpart = ipart + 1; nn < d_aa.M[ipart + 1]; nn++, jpart++)
      {
        RealType rV, d_rV_dr, d2_rV_dr2;
        rV         = rVsforce->splint(d_aa.r(nn), d_rV_dr, d2_rV_dr2);
        RealType V = rV * d_aa.rinv(nn);
        esum += Zat[jpart] * d_aa.rinv(nn) * rV;
        PosType grad = Zat[jpart] * Zat[ipart] * (d_rV_dr - V) * d_aa.rinv(nn) * d_aa.rinv(nn) * d_aa.dr(nn);
        forces[ipart] += grad;
        forces[jpart] -= grad;
      }
      //Accumulate pair sums...species charge for atom i.
      SR += Zat[ipart] * esum;
    }
#endif
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
  const DistanceTableData& d_aa(P.getDistTable(d_aa_ID));
  mRealType SR = 0.0;
  if (d_aa.DTType == DT_SOA)
  {
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
  }
  else
  { //this will be removed
#ifndef ENABLE_SOA
    for (int ipart = 0; ipart < NumCenters; ipart++)
    {
      mRealType esum = 0.0;
      for (int nn = d_aa.M[ipart], jpart = ipart + 1; nn < d_aa.M[ipart + 1]; nn++, jpart++)
      {
        //if(d_aa.r(nn)>=myRcut) continue;
        //esum += Zat[jpart]*AA->evaluate(d_aa.r(nn),d_aa.rinv(nn));
        esum += Zat[jpart] * d_aa.rinv(nn) * rVs->splint(d_aa.r(nn));
      }
      //Accumulate pair sums...species charge for atom i.
      SR += Zat[ipart] * esum;
    }
#endif
  }
  return SR;
}

CoulombPBCAA::Return_t CoulombPBCAA::evalLR(ParticleSet& P)
{
  mRealType res = 0.0;
  const StructFact& PtclRhoK(*(P.SK));
  if (PtclRhoK.SuperCellEnum == SUPERCELL_SLAB)
  {
    const DistanceTableData& d_aa(P.getDistTable(d_aa_ID));
    if (d_aa.DTType == DT_SOA)
    {
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
                                PtclRhoK.KLists.kshell, PtclRhoK.eikr[iat], PtclRhoK.eikr[jat]);
#endif
        res += Zat[iat] * u;
      }
    }
    else
    {
      //distance table handles jat>iat
      for (int iat = 0; iat < NumCenters; ++iat)
      {
        mRealType u = 0;
#if !defined(USE_REAL_STRUCT_FACTOR)
        for (int nn = d_aa.M[iat], jat = iat + 1; nn < d_aa.M[iat + 1]; ++nn, ++jat)
          u += Zat[jat] *
              AA->evaluate_slab(d_aa.dr(nn)[slab_dir], PtclRhoK.KLists.kshell, PtclRhoK.eikr[iat], PtclRhoK.eikr[jat]);
#endif
        res += Zat[iat] * u;
      }
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
        mRealType temp = AA->evaluate(PtclRhoK.KLists.kshell, PtclRhoK.rhok_r[spec1], PtclRhoK.rhok_i[spec1],
                                      PtclRhoK.rhok_r[spec2], PtclRhoK.rhok_i[spec2]);
#else
        mRealType temp = AA->evaluate(PtclRhoK.KLists.kshell, PtclRhoK.rhok[spec1], PtclRhoK.rhok[spec2]);
#endif
        if (spec2 == spec1)
          temp *= 0.5;
        res += Z1 * Zspec[spec2] * temp;
      } //spec2
    }   //spec1
  }
  return res;
}


CoulombPBCAA::Return_t CoulombPBCAA::evalConsts_orig(bool report)
{
  //LRHandlerType::BreakupBasisType &Basis(AA->Basis);
  //const Vector<RealType> &coefs(AA->coefs);
  RealType Consts = 0.0; // constant term
  //v_l(r=0) including correction due to the non-periodic direction
  RealType vl_r0 = AA->evaluateLR_r0();
  for (int spec = 0; spec < NumSpecies; spec++)
  {
    RealType z = Zspec[spec];
    RealType n = NofSpecies[spec];
    Consts -= 0.5 * vl_r0 * z * z * n;
  }
  if (report)
    app_log() << "   PBCAA self-interaction term " << Consts << std::endl;
  //Compute Madelung constant: this is not correct for general cases
  MC0 = 0.0;
  for (int i = 0; i < AA->Fk.size(); i++)
    MC0 += AA->Fk[i];
  MC0 = 0.5 * (MC0 - vl_r0);
  //Neutraling background term
  RealType vs_k0 = AA->evaluateSR_k0(); //v_s(k=0)
  for (int speca = 0; speca < NumSpecies; speca++)
  {
    RealType za = Zspec[speca];
    RealType na = NofSpecies[speca];
    Consts -= 0.5 * vs_k0 * za * na * za * na;
    for (int specb = speca + 1; specb < NumSpecies; specb++)
    {
      RealType zb = Zspec[specb];
      int nb      = NofSpecies[specb];
      Consts -= vs_k0 * za * zb * na * nb;
    }
  }
  if (report)
    app_log() << "   PBCAA total constant " << Consts << std::endl;
  //app_log() << "   MC0 of PBCAA " << MC0 << std::endl;
  return Consts;
}


CoulombPBCAA::Return_t CoulombPBCAA::evalSR_old(ParticleSet& P)
{
  const auto& d_aa = P.getDistTable(d_aa_ID);
  RealType SR      = 0.0;
#ifndef ENABLE_SOA
  for (int ipart = 0; ipart < NumCenters; ipart++)
  {
    RealType esum = 0.0;
    for (int nn = d_aa.M[ipart], jpart = ipart + 1; nn < d_aa.M[ipart + 1]; nn++, jpart++)
    {
      //if(d_aa.r(nn)>=myRcut) continue;
      //esum += Zat[jpart]*AA->evaluate(d_aa.r(nn),d_aa.rinv(nn));
      esum += Zat[jpart] * d_aa.rinv(nn) * rVs->splint(d_aa.r(nn));
    }
    //Accumulate pair sums...species charge for atom i.
    SR += Zat[ipart] * esum;
  }
#endif
  return SR;
}

CoulombPBCAA::Return_t CoulombPBCAA::evalLR_old(ParticleSet& P)
{
  RealType LR = 0.0;
  const StructFact& PtclRhoK(*(P.SK));
  for (int spec1 = 0; spec1 < NumSpecies; spec1++)
  {
    RealType Z1 = Zspec[spec1];
    for (int spec2 = spec1; spec2 < NumSpecies; spec2++)
    {
      RealType Z2 = Zspec[spec2];
#if defined(USE_REAL_STRUCT_FACTOR)
      RealType temp = AA->evaluate(PtclRhoK.KLists.kshell, PtclRhoK.rhok_r[spec1], PtclRhoK.rhok_i[spec1],
                                   PtclRhoK.rhok_r[spec2], PtclRhoK.rhok_i[spec2]);
#else
      RealType temp = AA->evaluate(PtclRhoK.KLists.kshell, PtclRhoK.rhok[spec1], PtclRhoK.rhok[spec2]);
#endif
      if (spec2 == spec1)
        LR += 0.5 * Z1 * Z2 * temp;
      else
        LR += Z1 * Z2 * temp;
    } //spec2
  }   //spec1
  //LR*=0.5;
  return LR;
}

CoulombPBCAA::Return_t CoulombPBCAA::evalConsts_old(bool report)
{
  //LRHandlerType::BreakupBasisType &Basis(AA->Basis);
  //const Vector<RealType> &coefs(AA->coefs);
  RealType Consts = 0.0, V0 = 0.0;
  //for(int n=0; n<coefs.size(); n++)
  //  V0 += coefs[n]*Basis.h(n,0.0); //For charge q1=q2=1
  V0 = AA->evaluateLR_r0();
  for (int spec = 0; spec < NumSpecies; spec++)
  {
    RealType z = Zspec[spec];
    RealType n = NofSpecies[spec];
    Consts += -V0 * 0.5 * z * z * n;
  }
  //V0 = Basis.get_rc()*Basis.get_rc()*0.5;
  //for(int n=0; n<Basis.NumBasisElem(); n++)
  //  V0 -= coefs[n]*Basis.hintr2(n);
  //V0 *= 2.0*TWOPI/Basis.get_CellVolume(); //For charge q1=q2=1
  V0 = AA->evaluateSR_k0();
  for (int spec = 0; spec < NumSpecies; spec++)
  {
    RealType z = Zspec[spec];
    int n      = NofSpecies[spec];
    Consts += -V0 * z * z * 0.5 * n * n;
  }
  //If we have more than one species in this particleset then there is also a
  //single AB term that should be added to the last constant...
  //=-Na*Nb*V0*Za*Zb
  //This accounts for the partitioning of the neutralizing background...
  for (int speca = 0; speca < NumSpecies; speca++)
  {
    RealType za = Zspec[speca];
    int na      = NofSpecies[speca];
    for (int specb = speca + 1; specb < NumSpecies; specb++)
    {
      RealType zb = Zspec[specb];
      int nb      = NofSpecies[specb];
      Consts += -V0 * za * zb * na * nb;
    }
  }
  if (report)
    app_log() << "   Constant of PBCAA " << Consts << std::endl;
  return Consts;
}


OperatorBase* CoulombPBCAA::makeClone(ParticleSet& qp, TrialWaveFunction& psi)
{
  if (is_active)
    return new CoulombPBCAA(qp, is_active, ComputeForces);
  else
    return new CoulombPBCAA(*this); //nothing needs to be re-evaluated
}
} // namespace qmcplusplus
