//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
    
    
#include "QMCHamiltonians/StressPBCAA.h"
#include "Particle/DistanceTable.h"
#include "Particle/DistanceTableData.h"
#include "Utilities/ProgressReportEngine.h"
#include <numeric>

namespace qmcplusplus
{

StressPBCAA::StressPBCAA(ParticleSet& ref, bool active) :
  AA(0), myGrid(0), rVs(0), FirstTime(true), myConst(0.0), ForceBase(ref,ref), Ps(ref), is_active(active)
{
  ReportEngine PRE("StressPBCAA","StressPBCAA");
  //save source tag
  SourceID=ref.tag();
  //create a distance table: just to get the table name
  DistanceTableData *d_aa = DistanceTable::add(ref);
  PtclRefName=d_aa->Name;
  initBreakup(ref);
  prefix="S_"+PtclRefName;
  app_log() << "  Maximum K shell " << AA->MaxKshell << std::endl;
  app_log() << "  Number of k vectors " << AA->Fk.size() << std::endl;
  if(!is_active)
  {
    d_aa->evaluate(ref);
    update_source(ref);
   app_log()<<"Evaluating Stress SymTensor::Long Range\n"; 
    sLR=evalLR(ref);
   app_log()<<"Short Range...\n";
    sSR=evalSR(ref);
    stress=sLR+sSR+myConst;
    
    //RealType eL(0.0), eS(0.0);
    //if (computeForces)
    //{
    //  forces = 0.0;
    //  eS=evalSRwithForces(ref);
    //  // 1.3978248322
    //  eL=evalLRwithForces(ref);
    //  // 2.130267378
    //}
    //else
    //{
    //  eL=evalLR(ref);
    //  eS=evalSR(ref);
    //}
    //NewValue=Value = eL+eS+myConst;
    //app_log() << "  Fixed Coulomb potential for " << ref.getName();
    //app_log() << "\n    e-e Madelung Const. =" << MC0
    //          << "\n    Vtot     =" << Value << std::endl;
    
  }
  app_log() << "  Stress SymTensor components for  " << ref.getName();
  app_log() << "\n    e-e Madelung Const. =\n" << MC0
            << "\n    Stot     =\n" << stress 
            << "\n    S_SR     =\n" << sSR   
            << "\n    S_LR     =\n" << sLR
            << "\n    S_Const  =\n" << myConst<< std::endl;
}

StressPBCAA:: ~StressPBCAA() { }


void StressPBCAA::update_source(ParticleSet& s)
{
 /* if(s.tag() == SourceID || s.parent() == SourceID)
  {
    RealType eL(0.0), eS(0.0);
    if (ComputeForces)
    {
      forces = 0.0;
      eS=evalSRwithForces(s);
      eL=evalLRwithForces(s);
    }
    else
    {
      eL=evalLR(s);
      eS=evalSR(s);
    }
    NewValue=Value = eL+eS+myConst;
  }*/
}

void StressPBCAA::resetTargetParticleSet(ParticleSet& P)
{
  if(is_active)
  {
    PtclRefName=P.DistTables[0]->Name;
    AA->resetTargetParticleSet(P);
  }
}
/*

void StressPBCAA::checkout_particle_arrays(TraceManager& tm)
{
//  V_sample = tm.checkout_real<1>(myName,Ps);
 // if(!is_active)
 //   spevaluate(Ps);
}

void StressPBCAA::delete_particle_arrays()
{
 // delete V_sample;
}
*/

StressPBCAA::Return_t
StressPBCAA::evaluate(ParticleSet& P)
{
  if(is_active)
  {
  //  if(tracing_particle_quantities)
 //     Value = spevaluate(P);
      stress = evalLR(P)+evalSR(P)+myConst;
  }
 // return Value;
 return 0.0;
}

/*
StressPBCAA::Return_t
StressPBCAA::spevaluate(ParticleSet& P)
{
  RealType  Vsr = 0.0;
  RealType  Vlr = 0.0;
  RealType& Vc  = myConst;
 // Array<RealType,1>& V_samp = *V_sample;
 // V_samp = 0.0;
  {
    //SR
    const DistanceTableData &d_aa(*P.DistTables[0]);
    RealType pairpot; //energy for single pair
    RealType z;
    for(int ipart=0; ipart<NumCenters; ipart++)
    {
      z = .5*Zat[ipart];
      for(int nn=d_aa.M[ipart],jpart=ipart+1; nn<d_aa.M[ipart+1]; nn++,jpart++)
      {
        pairpot = z*Zat[jpart]*d_aa.rinv(nn)*rVs->splint(d_aa.r(nn));
  //      V_samp(ipart)+=pairpot;
  //      V_samp(jpart)+=pairpot;
        Vsr += pairpot;
      }
    }
    Vsr *= 2.0;
  }
  {
    //LR
    const StructFact& PtclRhoK(*(P.SK));
    if(PtclRhoK.SuperCellEnum==SUPERCELL_SLAB)
    {
      APP_ABORT("StressPBCAA::spevaluate single particle traces have not been implemented for slab geometry");
    }
    else
    {
      //jtk mark: needs optimizations for USE_REAL_STRUCT_FACTOR
      RealType v1; //single particle energy
      RealType z;
      for(int i=0; i<NumCenters; i++)
      {
        z  = .5*Zat[i];
        v1 = 0.0;
        for(int s=0; s<NumSpecies; ++s){
#if defined(USE_REAL_STRUCT_FACTOR)
          v1 += z*Zspec[s]*AA->evaluate(PtclRhoK.KLists.kshell,PtclRhoK.rhok_r[s],PtclRhoK.rhok_i[s],PtclRhoK.eikr_r[i],PtclRhoK.eikr_i[i]);
#else
          v1 += z*Zspec[s]*AA->evaluate(PtclRhoK.KLists.kshell,PtclRhoK.rhok[s],PtclRhoK.eikr[i]);
#endif
        }
    //    V_samp(i)+=v1;
        Vlr += v1;
      }
    }
  }
  for(int i=0; i<V_sample->size(); ++i)
    V_samp(i)+=V_const(i);
  Value = Vsr + Vlr + Vc;
#if defined(TRACE_CHECK)
  RealType Vlrnow  = evalLR(P);
  RealType Vsrnow  = evalSR(P);
  RealType Vcnow   = myConst;
  RealType Vnow    = Vlrnow+Vsrnow+Vcnow;
  RealType Vsum    = V_samp.sum();
  RealType Vcsum   = V_const.sum();
  RealType Vsrold  = evalSR_old(P);
  RealType Vlrold  = evalLR_old(P);
  RealType Vcold   = evalConsts_old(false);
  RealType Vcorig  = evalConsts_orig(false);
  if(std::abs(Vsum-Vnow)>TraceManager::trace_tol)
  {
    app_log()<<"accumtest: StressPBCAA::evaluate()"<< std::endl;
    app_log()<<"accumtest:   tot:"<< Vnow << std::endl;
    app_log()<<"accumtest:   sum:"<< Vsum  << std::endl;
    APP_ABORT("Trace check failed");
  }
  if(std::abs(Vcsum-Vcnow)>TraceManager::trace_tol)
  {
    app_log()<<"accumtest: StressPBCAA::evalConsts()"<< std::endl;
    app_log()<<"accumtest:   tot:"<< Vcnow << std::endl;
    app_log()<<"accumtest:   sum:"<< Vcsum  << std::endl;
    APP_ABORT("Trace check failed");
  }
  if(std::abs(Vsrold-Vsrnow)>TraceManager::trace_tol)
  {
    app_log()<<"versiontest: StressPBCAA::evalSR()"<< std::endl;
    app_log()<<"versiontest:    old:"<< Vsrold << std::endl;
    app_log()<<"versiontest:    mod:"<< Vsrnow << std::endl;
    APP_ABORT("Trace check failed");
  }
  if(std::abs(Vlrold-Vlrnow)>TraceManager::trace_tol)
  {
    app_log()<<"versiontest: StressPBCAA::evalLR()"<< std::endl;
    app_log()<<"versiontest:    old:"<< Vlrold << std::endl;
    app_log()<<"versiontest:    mod:"<< Vlrnow << std::endl;
    APP_ABORT("Trace check failed");
  }
  if(std::abs(Vcold-Vcorig)>TraceManager::trace_tol ||
      std::abs(Vcnow-Vcorig)>TraceManager::trace_tol )
  {
    app_log()<<"versiontest: StressPBCAA::evalConsts()"<< std::endl;
    app_log()<<"versiontest:    old:"<< Vcold << std::endl;
    app_log()<<"versiontest:   orig:"<< Vcorig << std::endl;
    app_log()<<"versiontest:    mod:"<< Vcnow << std::endl;
    APP_ABORT("Trace check failed");
  }
#endif
  return Value;
}
*/
void StressPBCAA::initBreakup(ParticleSet& P)
{
  //SpeciesSet& tspecies(PtclRef->getSpeciesSet());
  SpeciesSet& tspecies(P.getSpeciesSet());
  //Things that don't change with lattice are done here instead of InitBreakup()
  ChargeAttribIndx = tspecies.addAttribute("charge");
  MemberAttribIndx = tspecies.addAttribute("membersize");
  NumCenters = P.getTotalNum();
  NumSpecies = tspecies.TotalNum;
 // V_const.resize(NumCenters);
  Zat.resize(NumCenters);
  Zspec.resize(NumSpecies);
  NofSpecies.resize(NumSpecies);
  for(int spec=0; spec<NumSpecies; spec++)
  {
    Zspec[spec] = tspecies(ChargeAttribIndx,spec);
    NofSpecies[spec] = static_cast<int>(tspecies(MemberAttribIndx,spec));
  }
  SpeciesID.resize(NumCenters);
  for(int iat=0; iat<NumCenters; iat++)
  {
    SpeciesID[iat]=P.GroupID[iat];
    Zat[iat] = Zspec[P.GroupID[iat]];
  }
  AA = LRCoulombSingleton::getDerivHandler(P);
  //AA->initBreakup(*PtclRef);
  myConst=evalConsts();
  myRcut=AA->get_rc();//Basis.get_rc();
  if(rVs==0)
  {
    rVs = LRCoulombSingleton::createSpline4RbyVs(AA,myRcut,myGrid);
  }
}


/*StressPBCAA::Return_t
StressPBCAA::evalLRwithForces(ParticleSet& P)
{
  RealType LR=0.0;
  const StructFact& PtclRhoK(*(P.SK));
  std::vector<TinyVector<RealType,DIM> > grad(P.getTotalNum());
  for(int spec2=0; spec2<NumSpecies; spec2++)
  {
    RealType Z2 = Zspec[spec2];
    for (int iat=0; iat<grad.size(); iat++)
      grad[iat] = TinyVector<RealType,DIM>(0.0);
    AA->evaluateGrad(P, P, spec2, Zat, grad);
    for (int iat=0; iat<grad.size(); iat++)
      forces[iat] += Z2*grad[iat];
  } //spec2
  return evalLR(P);
}*/



/*StressPBCAA::Return_t
StressPBCAA::evalSRwithForces(ParticleSet& P)
{
  const DistanceTableData *d_aa = P.DistTables[0];
  RealType SR=0.0;
  for(int ipart=0; ipart<NumCenters; ipart++)
  {
    RealType esum = 0.0;
    for(int nn=d_aa->M[ipart],jpart=ipart+1; nn<d_aa->M[ipart+1]; nn++,jpart++)
    {
      RealType rV, d_rV_dr, d2_rV_dr2;
      rV = rVs->splint(d_aa->r(nn), d_rV_dr, d2_rV_dr2);
      RealType V = rV *d_aa->rinv(nn);
      esum += Zat[jpart]*d_aa->rinv(nn)*rV;
      PosType grad = Zat[jpart]*Zat[ipart]*
                     (d_rV_dr - V)*d_aa->rinv(nn)*d_aa->rinv(nn)*d_aa->dr(nn);
      forces[ipart] += grad;
      forces[jpart] -= grad;
    }
    //Accumulate pair sums...species charge for atom i.
    SR += Zat[ipart]*esum;
  }
  return SR;
}*/


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
SymTensor<StressPBCAA::RealType, OHMMS_DIM>
StressPBCAA::evalConsts(bool report)
{
  SymTensor<RealType, OHMMS_DIM> Consts=0.0; // constant term
  RealType v1; //single particle energy
 // V_const = 0.0;
  //v_l(r=0) including correction due to the non-periodic direction
  SymTensor<RealType, OHMMS_DIM> vl_r0 = AA->evaluateLR_r0_dstrain();
  app_log()<<"   PBCAA vl_r0 = \n"<<vl_r0<< std::endl;
  for(int ipart=0; ipart<NumCenters; ipart++)
  {
    Consts += -.5*Zat[ipart]*Zat[ipart]*vl_r0;
  }
  
  if(report)
    app_log() << "   PBCAA self-interaction term\n " << Consts << std::endl;
  //Compute Madelung constant: this is not correct for general cases
  MC0 = 0.0;

  const StructFact& PtclRhoK(*(Ps.SK));
  const std::vector<PosType> &kpts = PtclRhoK.KLists.kpts_cart;
  for(int ks=0; ks<AA->dFk_dstrain.size(); ks++)
    {
        MC0+=AA->dFk_dstrain[ks];
    }
  app_log()<<"   PBCAA MC0 = "<< MC0<< std::endl;
  MC0 = 0.5*(MC0 - vl_r0);
  //Neutraling background term
  SymTensor<RealType, OHMMS_DIM> vs_k0=AA->evaluateSR_k0_dstrain(); //v_s(k=0)
  app_log()<<"    PBCAA Background Term:\n"<<vs_k0<< std::endl;
  for(int ipart=0; ipart<NumCenters; ipart++)
  {
    v1 = 0.0;
    for(int spec=0; spec<NumSpecies; spec++)
      v1 += -.5*Zat[ipart]*NofSpecies[spec]*Zspec[spec];
 //   v1 *= -.5*Zat[ipart]*vs_k0;
 //   V_const(ipart) += v1;
    Consts += v1*vs_k0;
  }
  if(report)
    app_log() << "   PBCAA total constant \n" << Consts << std::endl;
  //app_log() << "   MC0 of PBCAA " << MC0 << std::endl;
  return Consts;
}





SymTensor<StressPBCAA::RealType, OHMMS_DIM>
StressPBCAA::evalSR(ParticleSet& P)
{
  const DistanceTableData &d_aa(*P.DistTables[0]);
  SymTensor<RealType, OHMMS_DIM> SR=0.0;
  for(int ipart=0; ipart<NumCenters; ipart++)
  {
    SymTensor<RealType, OHMMS_DIM> esum=0.0;
    for(int nn=d_aa.M[ipart],jpart=ipart+1; nn<d_aa.M[ipart+1]; nn++,jpart++)
    {
      esum += Zat[jpart]* AA->evaluateSR_dstrain(d_aa.dr(nn), d_aa.r(nn));
      //app_log()<<"ipart = "<<ipart<<" nn="<<nn<< " dr = "<< d_aa.dr(nn)<< std::endl;
      //app_log()<<"srDf(r) = "<<AA->srDf(d_aa.r(nn),1.0/d_aa.r(nn))<< " part stress = \n"<<Zat[jpart]* AA->evaluateSR_dstrain(d_aa.dr(nn), d_aa.r(nn))<< std::endl;
    }
    //Accumulate pair sums...species charge for atom i.
    SR += Zat[ipart]*esum;
    //app_log()<<" SRpartsum=\n"<<Zat[ipart]*esum<< std::endl<< std::endl;
  }
  return SR;
  
}

SymTensor<StressPBCAA::RealType, OHMMS_DIM>
StressPBCAA::evalLR(ParticleSet& P)
{

  SymTensor<RealType, OHMMS_DIM> res;
  const StructFact& PtclRhoK(*(P.SK));
  
  SymTensor<RealType, OHMMS_DIM> temp;
    for(int spec1=0; spec1<NumSpecies; spec1++)
    {
      RealType Z1 = Zspec[spec1];
      for(int spec2=spec1; spec2<NumSpecies; spec2++)
      {
		//temp=0;
 #if !defined(USE_REAL_STRUCT_FACTOR)
        SymTensor<RealType, OHMMS_DIM> temp = AA->evaluateStress(PtclRhoK.KLists.kshell, PtclRhoK.rhok[spec1], PtclRhoK.rhok[spec2]);
 #else      
        SymTensor<RealType, OHMMS_DIM> temp = AA->evaluateStress(PtclRhoK.KLists.kshell, PtclRhoK.rhok_r[spec1], PtclRhoK.rhok_i[spec1], PtclRhoK.rhok_r[spec2], PtclRhoK.rhok_i[spec2]);
 #endif       
       // app_log()<<"WEEEE "<<temp<< std::endl;ls
       
        if(spec2==spec1)
          temp*=0.5;
        res += Z1*Zspec[spec2]*temp;
      } //spec2
    }//spec1
  
  return res;
}




QMCHamiltonianBase* StressPBCAA::makeClone(ParticleSet& qp, TrialWaveFunction& psi)
{
  if(is_active)
    return new StressPBCAA(qp, is_active);
  else
    return new StressPBCAA(*this);//nothing needs to be re-evaluated
}
}


