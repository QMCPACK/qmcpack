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
    
    
#include "QMCHamiltonians/ForceChiesaPBCAA.h"
#include "Particle/DistanceTable.h"
#include "Particle/DistanceTableData.h"
#include "Message/Communicate.h"
#include "Utilities/ProgressReportEngine.h"
#include "Numerics/DeterminantOperators.h"
#include "Numerics/MatrixOperators.h"
#include "OhmmsData/ParameterSet.h"
#include "OhmmsData/AttributeSet.h"

namespace qmcplusplus
{

ForceChiesaPBCAA::ForceChiesaPBCAA(ParticleSet& ions, ParticleSet& elns, bool firsttime):
  ForceBase(ions, elns), PtclA(ions), first_time(firsttime)
{
  ReportEngine PRE("ForceChiesaPBCAA","ForceChiesaPBCAA");
  myName = "Chiesa_Force_Base_PBCAB";
  prefix="FChiesaPBC";
  //Defaults for the chiesa S-wave polynomial filtering.  
  Rcut = 0.4;
  m_exp = 2;
  N_basis = 4;
  forces = 0.0;
  forces_ShortRange.resize(Nnuc);
  forces_ShortRange = 0.0;
  forces_IonIon=0.0;
  ions.turnOnPerParticleSK();
  //This sets up the long range breakups. 
  kcdifferent=false;
  myTableIndex=elns.addTable(ions,DT_SOA_PREFERRED);
  initBreakup(elns);
 // app_log()<< "IonIon Force" << std::endl;
 // app_log()<<forces_IonIon<< std::endl; 
  if (first_time==true)
  { 
   evaluateLR_AA();
 // app_log()<< "IonIon Force" << std::endl;
 // app_log()<<forces_IonIon<< std::endl; 
   evaluateSR_AA();
   app_log()<< "IonIon Force" << std::endl;
   app_log()<<forces_IonIon<< std::endl;
   first_time=false;
  }
 // forces=0.0;
 // evaluateLR(elns);
 // app_log()<<"LR eI FORCE\n";
 // app_log()<<forces<< std::endl;

 // evaluateSR(elns);
 // app_log()<<"LR+SR eI FORCE\n";
 // app_log()<<forces<< std::endl;
  
 
 // app_log() << "  Maximum K shell " << AB->MaxKshell << std::endl;
 // app_log() << "  Number of k vectors " << AB->Fk.size() << std::endl;
  
  ///////////////////////////////////////////////////////////////
}

void ForceChiesaPBCAA::InitMatrix()
{
  Sinv.resize(N_basis, N_basis);
  h.resize(N_basis);
  c.resize(N_basis);
  for(int k=0; k<N_basis; k++)
  {
    h[k] = std::pow(Rcut, (k+2))/static_cast<RealType>(k+2);
    for(int j=0; j<N_basis; j++)
    {
      Sinv(k,j) = std::pow(Rcut, (m_exp+k+j+3))/static_cast<RealType>(m_exp+k+j+3);
    }
  }
  // in Numerics/DeterminantOperators.h
  invert_matrix(Sinv, false);
  // in Numerics/MatrixOperators.h
  MatrixOperators::product(Sinv, h.data(), c.data());
}

void ForceChiesaPBCAA::initBreakup(ParticleSet& P)
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
//    if(totQ>std::numeric_limits<RealType>::epsilon())
//    {
//      LOGMSG("PBCs not yet finished for non-neutral cells");
//      OHMMS::Controller->abort();
//    }
  ////Test if the box sizes are same (=> kcut same for fixed dimcut)
  kcdifferent = (std::abs(PtclA.Lattice.LR_kc - P.Lattice.LR_kc) > std::numeric_limits<RealType>::epsilon());
  minkc = std::min(PtclA.Lattice.LR_kc,P.Lattice.LR_kc);
  //AB->initBreakup(*PtclB);
  //initBreakup is called only once
  //AB = LRCoulombSingleton::getHandler(*PtclB);
  AB = LRCoulombSingleton::getDerivHandler(P);
//  myConst=evalConsts();
  myRcut=AB->get_rc();//Basis.get_rc();
  // create the spline function for the short-range part assuming pure potential
 // if(V0==0)
 // {
 //   V0 = LRCoulombSingleton::createSpline4RbyVs(AB,myRcut,myGrid);
 //   if(Vat.size())
 //   {
 //     app_log() << "  Vat is not empty. Something is wrong" << std::endl;
 //     OHMMS::Controller->abort();
 //   }
 //   Vat.resize(NptclA,V0);
 //   Vspec.resize(NumSpeciesA,0);//prepare for PP to overwrite it
 // }
}

void ForceChiesaPBCAA::evaluateLR(ParticleSet& P)
{
  const StructFact& RhoKA(*(PtclA.SK));
  const StructFact& RhoKB(*(P.SK));
  //app_log()<<"Calculate Long Range e-I forces"<< std::endl;
  std::vector<TinyVector<RealType,DIM> > grad(PtclA.getTotalNum());
  for(int j=0; j<NumSpeciesB; j++)
  {
    for (int iat=0; iat<grad.size(); iat++)
      grad[iat] = TinyVector<RealType,DIM>(0.0, 0.0, 0.0);
    AB->evaluateGrad(PtclA, P, j, Zat, grad);
    for (int iat=0; iat<grad.size(); iat++){
      forces[iat] += Qspec[j]*grad[iat];
    //  app_log()<<"Qspec["<<j<<"] = "<<Qspec[j]<< std::endl;
  //    app_log()<<"Grad["<<iat<<"] = "<<grad[iat]<< std::endl;
    }
  } // electron species
  

}

void ForceChiesaPBCAA::evaluateSR(ParticleSet& P)
{
  const DistanceTableData &d_ab(*P.DistTables[myTableIndex]);
  if(d_ab.DTType == DT_SOA)
  {
    for(size_t jat=0; jat<NptclB; ++jat)
    {
      const RealType* restrict dist=d_ab.Distances[jat];

      for(size_t iat=0; iat<NptclA; ++iat)
      {
        const RealType r = dist[iat];
        const RealType rinv = RealType(1)/r;
        RealType g_f = g_filter(r);
        RealType V = -AB->srDf(r,rinv);
        PosType drhat = rinv*d_ab.Displacements[jat][iat];
        forces[iat] += g_f*Zat[iat]*Qat[jat]*V*drhat;
      }
    }
    //for(size_t a=0; a<NptclA; ++a)
    //  std::cout << "debug evaluateSR " << forces[a] << std::endl;
  }
  else
  {
    //RealType res=0.0;
    //Loop over distinct eln-ion pairs
    for(int iat=0; iat<NptclA; iat++)
    {
      //RealType esum = 0.0;
      //app_log()<<"Long Range force calculation for ion..."<< std::endl;
      for(int nn=d_ab.M[iat], jat=0; nn<d_ab.M[iat+1]; ++nn,++jat)
      {
        RealType V;
        RealType g_f=g_filter(d_ab.r(nn));
        //rV = rVs->splint(d_ab.r(nn), d_rV_dr, d2_rV_dr2);
        V = -AB->srDf(d_ab.r(nn),d_ab.rinv(nn));
       // std::stringstream wee;
      //  wee<<"srDf() #"<<omp_get_thread_num()<<" V= "<<V<<" "<<iat<<" "<<nn<< std::endl;
      //  std::cout <<wee.str();
        PosType drhat = d_ab.rinv(nn) * d_ab.dr(nn);
        //esum += Qat[jat]*d_ab.rinv(nn)*rV;
      //  app_log()<<"iat="<<iat<<" elec="<<jat<< std::endl;
      //  app_log()<<"dr = "<<d_ab.dr(nn)<< std::endl;
      //  app_log()<<"force = "<<-g_f*Zat[iat]*Qat[jat] * V * drhat<< std::endl;
        forces[iat] += -g_f*Zat[iat]*Qat[jat] * V * drhat;
      }
      //std::cout << "debug evaluateSR " << forces[iat] << std::endl;
    }
  }
}

void ForceChiesaPBCAA::evaluateSR_AA()
{
  const DistanceTableData &d_aa(*PtclA.DistTables[0]);
  if(d_aa.DTType == DT_SOA)
  {
    for(size_t ipart=1; ipart<NptclA; ipart++)
    {
      const RealType* restrict dist=d_aa.Distances[ipart];
      for(size_t jpart=0; jpart<ipart; ++jpart)
      {
        RealType V = -AB->srDf(dist[jpart],RealType(1)/dist[jpart]);
        PosType grad = -Zat[jpart]*Zat[ipart]*V*dist[jpart]*d_aa.Displacements[ipart][jpart];
        forces_IonIon[ipart] += grad;
        forces_IonIon[jpart] -= grad;
      }
    }
    for(size_t ipart=0; ipart<NptclA; ++ipart)
      std::cout << "debug evaluateSR_AA " << forces_IonIon[ipart] << std::endl;
  }
  else
  {
    //RealType res=0.0;
    //Loop over distinct eln-ion pairs
    for(int ipart=0; ipart<NptclA; ipart++)
    {
      for(int nn=d_aa.M[ipart],jpart=ipart+1; nn<d_aa.M[ipart+1]; nn++,jpart++)
      {
        RealType V = -AB->srDf(d_aa.r(nn),d_aa.rinv(nn));
        PosType grad = -Zat[jpart]*Zat[ipart]*V*d_aa.rinv(nn)*d_aa.dr(nn);
        forces_IonIon[ipart] += grad;
        forces_IonIon[jpart] -= grad;
     //   app_log()<<"ShortRange Ion Ion component"<< std::endl;
     //   app_log() <<"grad[" <<ipart<< "] = "<<grad<< std::endl;
     //   app_log() <<"Zat[" <<ipart<< "] = "<<Zat[ipart]<< std::endl;
      }
      std::cout << "debug evaluateSR_AA " << forces_IonIon[ipart] << std::endl;
    }
  }
}

void ForceChiesaPBCAA::evaluateLR_AA()
{
  const StructFact& PtclRhoK(*(PtclA.SK));
  std::vector<TinyVector<RealType,DIM> > grad(PtclA.getTotalNum());
  for(int spec2=0; spec2<NumSpeciesA; spec2++)
  {
    RealType Z2 = Zspec[spec2];
    for (int iat=0; iat<grad.size(); iat++)
      grad[iat] = TinyVector<RealType,DIM>(0.0);
    AB->evaluateGrad(PtclA, PtclA, spec2, Zat, grad);

    for (int iat=0; iat<grad.size(); iat++){
     // app_log()<<"Long Range Ion Ion Component."<< std::endl;
     // app_log() <<"grad[" <<iat<< "] = "<<grad[iat]<< std::endl;
     // app_log() <<"Zat[" <<iat<< "] = "<<Zat[iat]<< std::endl;
      forces_IonIon[iat] += Z2*grad[iat];
      
  }
  } //spec2
  
  
}


ForceChiesaPBCAA::Return_t
ForceChiesaPBCAA::evaluate(ParticleSet& P)
{
  //forces = forces_IonIon;
  forces=0.0;
  evaluateLR(P);
//  app_log()<<"LR eI FORCE\n";
//  app_log()<<forces<< std::endl;

  evaluateSR(P);
//  app_log()<<"LR+SR eI FORCE\n";
//  app_log()<<forces<< std::endl;
  if(addionion==true) forces=forces+forces_IonIon;  
  return 0.0;
}

ForceChiesaPBCAA::Return_t ForceChiesaPBCAA::g_filter(RealType r)
{
	if(r>=Rcut)
    {
       return 1.0;
    }
    else
    {
       RealType g_q=0.0;
       for (int q=0; q<N_basis; q++)
       {
	      g_q += c[q]*std::pow(r,m_exp+q+1);
	   }    
       
       return g_q;
    }
}

bool ForceChiesaPBCAA::put(xmlNodePtr cur)
{
  std::string ionionforce("yes");
  OhmmsAttributeSet attr;
  attr.add(prefix, "name");
  attr.add(ionionforce, "addionion");
  attr.put(cur);
  addionion = (ionionforce=="yes") || (ionionforce == "true");
  app_log() << "ionionforce = "<<ionionforce<< std::endl;
  app_log() << "addionion="<<addionion<< std::endl;
  app_log() << "FirstTime= "<<FirstTime<< std::endl;
  ParameterSet fcep_param_set;
  fcep_param_set.add(Rcut, "rcut","real");
  fcep_param_set.add(N_basis, "nbasis", "int");
  fcep_param_set.add(m_exp, "weight_exp", "int");
  fcep_param_set.put(cur);
  app_log() <<"    ForceChiesaPBCAA Parameters"<< std::endl;
  app_log() <<"        ForceChiesaPBCAA::Rcut="<<Rcut<< std::endl;
  app_log() <<"        ForceChiesaPBCAA::N_basis="<<N_basis<< std::endl;
  app_log() <<"        ForceChiesaPBCAA::m_exp="<<m_exp<< std::endl;
  InitMatrix();
  return true;
}

void ForceChiesaPBCAA::resetTargetParticleSet(ParticleSet& P)
{
  int tid=P.addTable(PtclA,DT_AOS);
  if(tid != myTableIndex)
  {
    APP_ABORT("ForceChiesaPBCAA::resetTargetParticleSet found inconsistent table index");
  }
  AB->resetTargetParticleSet(P);
}

void ForceChiesaPBCAA::addObservables(PropertySetType& plist, BufferType& collectables)
{
  myIndex=plist.add(myName.c_str());
//  if (ComputeForces)
    addObservablesF(plist);
}

QMCHamiltonianBase* ForceChiesaPBCAA::makeClone(ParticleSet& qp, TrialWaveFunction& psi)
{
//  ForceChiesaPBCAA* tmp = new ForceChiesaPBCAA(*this);
  ForceChiesaPBCAA* tmp= new ForceChiesaPBCAA(PtclA,qp,false);
  tmp->Rcut=Rcut; // parameter: radial distance within which estimator is used
  tmp->m_exp=m_exp; // parameter: exponent in polynomial fit
  tmp->N_basis=N_basis; // parameter: size of polynomial basis set
  tmp->Sinv.resize(N_basis,N_basis);
  tmp->Sinv=Sinv; // terms in fitting polynomial
  tmp->h.resize(N_basis);
  tmp->h=h; // terms in fitting polynomial
  tmp->c.resize(N_basis);
  tmp->c=c; // polynomial coefficients
  tmp->addionion=addionion;
  tmp->initBreakup(qp);

  return tmp;
}
}

//  void ForceChiesaPBCAA::addObservables(PropertySetType& plist) {
//    //cerr << "ForceBase::addObs sound off" << std::endl;
//    //obsName << myName << "0_x";
//    //myIndex = plist.add(obsName.str());
//    //obsName.clear();
//    mySize = Nnuc*OHMMS_DIM;
//    //cerr << "ForceBase mySize is " << Nnuc << " * " << OHMMS_DIM << " = " << mySize << std::endl;
//    checkInit = true;
//    if(myIndex<0) myIndex=plist.size();
//    int tmpIndex;
//    bool firstTime = true;
//    for(int iat=0; iat<Nnuc; iat++) {
//      for(int x=0; x<OHMMS_DIM; x++) {
//        std::ostringstream obsName;
//        obsName << "HFCep_" << iat << "_" << x;
//        tmpIndex = plist.add(obsName.str());
//        //if(firstTime) {
//        //  firstTime = false;
//        //  myIndex = tmpIndex;// + 3;
//        //}
//        std::cerr << iat << ", " << x << " stored at " << tmpIndex << std::endl;
//      }
//    }
//    std::cerr << "AddObs myIndex is " << myIndex << " last " << tmpIndex << std::endl;
//  }
//  // debugging version only z component
//  //void ForceChiesaPBCAA::addObservables(PropertySetType& plist) {
//  //  mySize = Nnuc*OHMMS_DIM;
//  //  checkInit = true;
//  //  int tmpIndex;
//  //  bool firstTime = true;
//  //  for(int iat=0; iat<Nnuc; iat++) {
//  //    //for(int x=0; x<OHMMS_DIM; x++) {
//  //      std::ostringstream obsName;
//  //      obsName << "HFCep_" << iat << "_Z_sr";
//  //      tmpIndex = plist.add(obsName.str());
//  //      if(firstTime) {
//  //        firstTime = false;
//  //        myIndex = tmpIndex;
//  //        std::cerr << "ForceChiesaPBCAA addObs setting myindex " << myIndex << std::endl;
//  //      }
//  //      std::ostringstream obsName2;
//  //      obsName2 << "HFCep_" << iat << "_Z_lr";
//  //      tmpIndex = plist.add(obsName2.str());
//  //      //cerr << iat << ", " << x << " stored at " << tmpIndex << std::endl;
//  //    //}
//  //  }
//  //}
//
//  //// overriding base class definition to print out short and long range contribs
//  //void ForceChiesaPBCAA::setObservables(PropertySetType& plist) {
//  //  //cerr << "ForceBase::setObs storing forces ";
//  //  int index = myIndex;
//  //  for(int iat=0; iat<Nnuc; iat++) {
//  //    // HACK only z component
//  //    //for(int x=0; x<OHMMS_DIM; x++) {
//  //    plist[index] = forces(iat,2);
//  //    //cerr << index << ": " << plist[index] << "; ";
//  //    index++;
//  //    plist[index] = storeForces(iat,2);
//  //    //cerr << index << ": " << plist[index] << "; ";
//  //    index++;
//  //    //}
//  //  }
//  //  //cerr << std::endl;
//  //}
//  void ForceChiesaPBCAA::setObservables(PropertySetType& plist) {
//    ///cerr << "ForceBase::setObs storing forces";
//    int index = myIndex;
//    for(int iat=0; iat<Nnuc; iat++) {
//      for(int x=0; x<OHMMS_DIM; x++) {
//        plist[index] = forces(iat,x) + storeForces(iat,x);
//        //cerr << " " << index << ":" << plist[index];
//        index++;
//      }
//    }
//    //cerr << std::endl;
//  }

