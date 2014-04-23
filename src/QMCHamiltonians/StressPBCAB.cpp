//////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Jeongnim Kim and Kris Delaney
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#include "QMCHamiltonians/StressPBCAB.h"
#include "Particle/DistanceTable.h"
#include "Particle/DistanceTableData.h"
#include "Message/Communicate.h"
#include "Utilities/ProgressReportEngine.h"

namespace qmcplusplus
{

StressPBCAB::StressPBCAB(ParticleSet& ions, ParticleSet& elns,
                           bool computeForces):
  PtclA(ions), myConst(0.0), myGrid(0),V0(0),ComputeForces(computeForces),
  ForceBase (ions, elns), MaxGridPoints(10000),Pion(ions),Peln(elns)
{
  // if (ComputeForces)
  // 	InitVarReduction (0.5, 0, 3);
  ReportEngine PRE("StressPBCAB","StressPBCAB");
  //Use singleton pattern
  //AB = new LRHandlerType(ions);
  myTableIndex=elns.addTable(ions);
  initBreakup(elns);
  prefix="S_AB";
  app_log() << "  Maximum K shell " << AB->MaxKshell << endl;
  app_log() << "  Number of k vectors " << AB->Fk.size() << endl;
  is_active=true;
}

QMCHamiltonianBase* StressPBCAB::makeClone(ParticleSet& qp, TrialWaveFunction& psi)
{
  StressPBCAB* myclone=new StressPBCAB(PtclA,qp,ComputeForces);
  myclone->FirstForceIndex = FirstForceIndex;
  if(myGrid)
    myclone->myGrid=new GridType(*myGrid);
  for(int ig=0; ig<Vspec.size(); ++ig)
  {
    if(Vspec[ig])
    {
      RadFunctorType* apot=Vspec[ig]->makeClone();
      myclone->Vspec[ig]=apot;
      for(int iat=0; iat<PtclA.getTotalNum(); ++iat)
      {
        if(PtclA.GroupID[iat]==ig)
          myclone->Vat[iat]=apot;
      }
    }
  }
  return myclone;
}

StressPBCAB:: ~StressPBCAB()
{
  //probably need to clean up
}

void StressPBCAB::resetTargetParticleSet(ParticleSet& P)
{
  int tid=P.addTable(PtclA);
  if(tid != myTableIndex)
  {
    APP_ABORT("StressPBCAB::resetTargetParticleSet found inconsistent table index");
  }
  AB->resetTargetParticleSet(P);
}


StressPBCAB::Return_t
StressPBCAB::evaluate(ParticleSet& P)
{
  if (is_active)
  {
    //forces = 0.0;
    stress = evalLR(P) + evalSR(P) +myConst;
  }


  return 0.0;
}






/** Evaluate the background term. Other constants are handled by AA potentials.
 *
 * \f$V_{bg}^{AB}=-\sum_{\alpha}\sum_{\beta} N^{\alpha} N^{\beta} q^{\alpha} q^{\beta} v_s(k=0) \f$
 * @todo Here is where the charge system has to be handled.
 */
SymTensor<StressPBCAB::RealType, OHMMS_DIM>
StressPBCAB::evalConsts(bool report)
{
  int nelns = Peln.getTotalNum();
  int nions = Pion.getTotalNum();
 // Ve_const.resize(nelns);
 // Vi_const.resize(nions);
 // Ve_const = 0.0;
 // Vi_const = 0.0;
  SymTensor<RealType, OHMMS_DIM> tmpconsts=0.0;
  SymTensor<RealType, OHMMS_DIM> vs_k0 = AB->evaluateSR_k0_dstrain();
  RealType v1;
 // SymTensor<RealType, OHMMS_DIM> v1; //single particle energy
  for(int i=0; i<nelns; ++i)
  {
    v1=0.0;
    for(int s=0; s<NumSpeciesA; s++)
      v1 += NofSpeciesA[s]*Zspec[s];
   // v1 *= -.5*Qat[i]*vs_k0;
    //Ve_const(i) = v1;
    tmpconsts += -.5*Qat[i]*vs_k0*v1;
  }
  for(int i=0; i<nions; ++i)
  {
    v1=0.0;
    for(int s=0; s<NumSpeciesB; s++)
      v1 += NofSpeciesB[s]*Qspec[s];
   // v1 *= -.5*Zat[i]*vs_k0;
   // Vi_const(i) = v1;
    tmpconsts += -.5*Zat[i]*vs_k0*v1;
  }
  //if(report)
    app_log() << "   Constant of PBCAB " << tmpconsts << endl;
  return tmpconsts;
} 






SymTensor<StressPBCAB::RealType, OHMMS_DIM>
StressPBCAB::evalSR(ParticleSet& P)
{
  const DistanceTableData &d_ab(*P.DistTables[myTableIndex]);
  SymTensor<RealType, OHMMS_DIM> res=0.0;
  //Loop over distinct eln-ion pairs
  for(int iat=0; iat<NptclA; iat++)
  {
    SymTensor<RealType, OHMMS_DIM> esum = 0.0;
   // RadFunctorType* rVs=Vat[iat];
    for(int nn=d_ab.M[iat], jat=0; nn<d_ab.M[iat+1]; ++nn,++jat)
    {
      // if(d_ab.r(nn)>=(myRcut-0.1)) continue;
      esum += Qat[jat]*AB->evaluateSR_dstrain(d_ab.dr(nn), d_ab.r(nn));
    }
    //Accumulate pair sums...species charge for atom i.
    res += Zat[iat]*esum;
  }
   app_log()<<"\nEvaluateSR_AA()_working = \n"<<res<<endl<<endl;
  return res;
}


SymTensor<StressPBCAB::RealType, OHMMS_DIM>
StressPBCAB::evalLR(ParticleSet& P)
{
 // const int slab_dir=OHMMS_DIM-1;
  SymTensor<RealType, OHMMS_DIM> res=0.0;
  const StructFact& RhoKA(*(PtclA.SK));
  const StructFact& RhoKB(*(P.SK));
//  if(RhoKA.SuperCellEnum==SUPERCELL_SLAB)
//  {
 //   const DistanceTableData &d_ab(*P.DistTables[myTableIndex]);
//    for(int iat=0; iat<NptclA; ++iat)
//    {
//      RealType u=0;
//#if !defined(USE_REAL_STRUCT_FACTOR)
//      for(int nn=d_ab.M[iat], jat=0; nn<d_ab.M[iat+1]; ++nn,++jat)
//        u += Qat[jat]*AB->evaluate_slab(d_ab.dr(nn)[slab_dir], RhoKA.KLists.kshell, RhoKA.eikr[iat], RhoKB.eikr[jat]);
//#endif
//      res += Zat[iat]*u;
//    }
//  }
//  else
//  {
    for(int i=0; i<NumSpeciesA; i++)
    {
      SymTensor<RealType, OHMMS_DIM> esum;
      esum=0.0;
      for(int j=0; j<NumSpeciesB; j++)
      {
#if defined(USE_REAL_STRUCT_FACTOR)
        esum += Qspec[j]*AB->evaluateStress(RhoKA.KLists.kshell
                                      , RhoKA.rhok_r[i],RhoKA.rhok_i[i] , RhoKB.rhok_r[j],RhoKB.rhok_i[j]);
#else
        esum += Qspec[j]*AB->evaluateStress(RhoKA.KLists.kshell, RhoKA.rhok[i],RhoKB.rhok[j]);
       // AA->evaluateStress(PtclRhoK.KLists.kshell, PtclRhoK.rhok[spec1], PtclRhoK.rhok[spec2], temp);
#endif
      } //speceln
      app_log()<<"\n   esum stressPBCAB = \n"<<esum<<endl<<endl;
      res += Zspec[i]*esum;
    }
//  }//specion
  app_log()<<"\nEvaluateLR_AB()_working = \n"<<res<<endl<<endl;
  return res;
}



void StressPBCAB::initBreakup(ParticleSet& P)
{
  SpeciesSet& tspeciesA(PtclA.getSpeciesSet());
  SpeciesSet& tspeciesB(P.getSpeciesSet());
  int ChargeAttribIndxA = tspeciesA.addAttribute("charge");
  int MemberAttribIndxA = tspeciesA.addAttribute("membersize");
  int ChargeAttribIndxB = tspeciesB.addAttribute("charge");
  int MemberAttribIndxB = tspeciesB.addAttribute("membersize");
  NptclA = PtclA.getTotalNum();
  NptclB = P.getTotalNum();
  NumSpeciesA = tspeciesA.TotalNum;
  NumSpeciesB = tspeciesB.TotalNum;
  //Store information about charges and number of each species
  Zat.resize(NptclA);
  Zspec.resize(NumSpeciesA);
  Qat.resize(NptclB);
  Qspec.resize(NumSpeciesB);
  NofSpeciesA.resize(NumSpeciesA);
  NofSpeciesB.resize(NumSpeciesB);
  for(int spec=0; spec<NumSpeciesA; spec++)
  {
    Zspec[spec] = tspeciesA(ChargeAttribIndxA,spec);
    NofSpeciesA[spec] = static_cast<int>(tspeciesA(MemberAttribIndxA,spec));
  }
  for(int spec=0; spec<NumSpeciesB; spec++)
  {
    Qspec[spec] = tspeciesB(ChargeAttribIndxB,spec);
    NofSpeciesB[spec] = static_cast<int>(tspeciesB(MemberAttribIndxB,spec));
  }
  RealType totQ=0.0;
  for(int iat=0; iat<NptclA; iat++)
    totQ+=Zat[iat] = Zspec[PtclA.GroupID[iat]];
  for(int iat=0; iat<NptclB; iat++)
    totQ+=Qat[iat] = Qspec[P.GroupID[iat]];
    
    
  for (int i=0; i<Zspec.size(); i++) app_log()<<"Zspec["<<i<<"]="<<Zspec[i]<<endl;
  for (int i=0; i<Qspec.size(); i++) app_log()<<"Qspec["<<i<<"]="<<Qspec[i]<<endl;
  for(int iat=0; iat<NptclA; iat++)
		app_log()<<"Zat["<<iat<<"]="<<Zat[iat]<<endl;
		
  for(int iat=0; iat<NptclB; iat++)
		app_log()<<"Qat["<<iat<<"]="<<Qat[iat]<<endl;      
//    if(totQ>numeric_limits<RealType>::epsilon())
//    {
//      LOGMSG("PBCs not yet finished for non-neutral cells");
//      OHMMS::Controller->abort();
//    }
  ////Test if the box sizes are same (=> kcut same for fixed dimcut)
  kcdifferent = (std::abs(PtclA.Lattice.LR_kc - P.Lattice.LR_kc) > numeric_limits<RealType>::epsilon());
  minkc = std::min(PtclA.Lattice.LR_kc,P.Lattice.LR_kc);
  //AB->initBreakup(*PtclB);
  //initBreakup is called only once
  //AB = LRCoulombSingleton::getHandler(*PtclB);
  AB = LRCoulombSingleton::getDerivHandler(P);
  myConst=evalConsts();
  myRcut=AB->get_rc();//Basis.get_rc();
  // create the spline function for the short-range part assuming pure potential
  if(V0==0)
  {
    V0 = LRCoulombSingleton::createSpline4RbyVs(AB,myRcut,myGrid);
    if(Vat.size())
    {
      app_log() << "  Vat is not empty. Something is wrong" << endl;
      OHMMS::Controller->abort();
    }
    Vat.resize(NptclA,V0);
    Vspec.resize(NumSpeciesA,0);//prepare for PP to overwrite it
  }
}


/*
StressPBCAB::Return_t
StressPBCAB::registerData(ParticleSet& P, BufferType& buffer)
{
  P.SK->DoUpdate=true;
  SRpart.resize(NptclB);
  LRpart.resize(NptclB);
  Value=evaluateForPyP(P);
  buffer.add(SRpart.begin(),SRpart.end());
  buffer.add(LRpart.begin(),LRpart.end());
  buffer.add(Value);
  return Value;
}


StressPBCAB::Return_t
StressPBCAB::updateBuffer(ParticleSet& P, BufferType& buffer)
{
  Value=evaluateForPyP(P);
  buffer.put(SRpart.begin(),SRpart.end());
  buffer.put(LRpart.begin(),LRpart.end());
  buffer.put(Value);
  return Value;
}

void StressPBCAB::copyFromBuffer(ParticleSet& P, BufferType& buffer)
{
  buffer.get(SRpart.begin(),SRpart.end());
  buffer.get(LRpart.begin(),LRpart.end());
  buffer.get(Value);
}

void StressPBCAB::copyToBuffer(ParticleSet& P, BufferType& buffer)
{
  buffer.put(SRpart.begin(),SRpart.end());
  buffer.put(LRpart.begin(),LRpart.end());
  buffer.put(Value);
}*/
/*
StressPBCAB::Return_t
StressPBCAB::evaluateForPyP(ParticleSet& P)
{
  Return_t res=myConst;
#if defined(USE_REAL_STRUCT_FACTOR)
  APP_ABORT("StressPBCAB::evaluateForPyP(ParticleSet& P)");
#else
  SRpart=0.0;
  const DistanceTableData* d_ab=P.DistTables[myTableIndex];
  for(int iat=0; iat<NptclA; ++iat)
  {
    RealType z=Zat[iat];
    RadFunctorType* rVs=Vat[iat];
    for(int nn=d_ab->M[iat], jat=0; nn<d_ab->M[iat+1]; nn++,jat++)
    {
      RealType e=z*Qat[jat]*d_ab->rinv(nn)*rVs->splint(d_ab->r(nn));
      SRpart[jat]+=e;
      res+=e;
    }
  }
  LRpart=0.0;
  const StructFact& RhoKA(*(PtclA.SK));
  const StructFact& RhoKB(*(P.SK));
  // const StructFact& RhoKB(*(PtclB->SK));
  for(int i=0; i<NumSpeciesA; i++)
  {
    RealType z=Zspec[i];
    for(int jat=0; jat<P.getTotalNum(); ++jat)
    {
      RealType e=z*Qat[jat]*AB->evaluate(RhoKA.KLists.kshell, RhoKA.rhok[i],RhoKB.eikr[jat]);
      LRpart[jat]+=e;
      res+=e;
    }
  }
#endif
  return res;
}


StressPBCAB::Return_t
StressPBCAB::evaluatePbyP(ParticleSet& P, int active)
{
#if defined(USE_REAL_STRUCT_FACTOR)
  APP_ABORT("StressPBCAB::evaluatePbyP(ParticleSet& P, int active)");
#else
  const std::vector<DistanceTableData::TempDistType> &temp(P.DistTables[myTableIndex]->Temp);
  RealType q=Qat[active];
  SRtmp=0.0;
  for(int iat=0; iat<NptclA; ++iat)
  {
    SRtmp+=Zat[iat]*q*temp[iat].rinv1*Vat[iat]->splint(temp[iat].r1);
  }
  LRtmp=0.0;
  const StructFact& RhoKA(*(PtclA.SK));
  //const StructFact& RhoKB(*(PtclB->SK));
  const StructFact& RhoKB(*(P.SK));
  for(int i=0; i<NumSpeciesA; i++)
    LRtmp+=Zspec[i]*q*AB->evaluate(RhoKA.KLists.kshell, RhoKA.rhok[i],RhoKB.eikr_temp.data());
#endif
  return NewValue=Value+(SRtmp-SRpart[active])+(LRtmp-LRpart[active]);
  //return NewValue=Value+(SRtmp-SRpart[active]);
}

void StressPBCAB::acceptMove(int active)
{
  SRpart[active]=SRtmp;
  LRpart[active]=LRtmp;
  Value=NewValue;
}


StressPBCAB::Return_t
StressPBCAB::evalLRwithForces(ParticleSet& P)
{
  const StructFact& RhoKA(*(PtclA.SK));
  const StructFact& RhoKB(*(P.SK));
  vector<TinyVector<RealType,DIM> > grad(PtclA.getTotalNum());
  for(int j=0; j<NumSpeciesB; j++)
  {
    for (int iat=0; iat<grad.size(); iat++)
      grad[iat] = TinyVector<RealType,DIM>(0.0, 0.0, 0.0);
    AB->evaluateGrad(PtclA, P, j, Zat, grad);
    for (int iat=0; iat<grad.size(); iat++)
      forces[iat] += Qspec[j]*grad[iat];
  } // electron species
  return evalLR(P);
}

StressPBCAB::Return_t
StressPBCAB::evalSRwithForces(ParticleSet& P)
{
  const DistanceTableData &d_ab(*P.DistTables[myTableIndex]);
  RealType res=0.0;
  //Loop over distinct eln-ion pairs
  for(int iat=0; iat<NptclA; iat++)
  {
    RealType esum = 0.0;
    RadFunctorType* rVs=Vat[iat];
    for(int nn=d_ab.M[iat], jat=0; nn<d_ab.M[iat+1]; ++nn,++jat)
    {
      RealType rV, d_rV_dr, d2_rV_dr2, V;
      rV = rVs->splint(d_ab.r(nn), d_rV_dr, d2_rV_dr2);
      V = rV *d_ab.rinv(nn);
      PosType drhat = d_ab.rinv(nn) * d_ab.dr(nn);
      esum += Qat[jat]*d_ab.rinv(nn)*rV;
      forces[iat] += Zat[iat]*Qat[jat] * //g(d_ab.r(nn)) *
                     (d_rV_dr - V)*d_ab.rinv(nn) *drhat;
    }
    //Accumulate pair sums...species charge for atom i.
    res += Zat[iat]*esum;
  }
  return res;
}*/
}

/***************************************************************************
 * $RCSfile$   $Author: jtkrogel $
 * $Revision: 5976 $   $Date: 2013-09-13 13:39:44 -0500 (Fri, 13 Sep 2013) $
 * $Id: StressPBCAB.cpp 5976 2013-09-13 18:39:44Z jtkrogel $
 ***************************************************************************/

