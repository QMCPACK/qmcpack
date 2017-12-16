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
    
    
#include "QMCHamiltonians/StressPBC.h"
#include "Particle/DistanceTable.h"
#include "Particle/DistanceTableData.h"
#include "Message/Communicate.h"
#include "Utilities/ProgressReportEngine.h"
#include "Numerics/DeterminantOperators.h"
#include "Numerics/MatrixOperators.h"
#include "OhmmsData/ParameterSet.h"
#include "OhmmsData/AttributeSet.h"
#include "QMCWaveFunctions/OrbitalBase.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"

namespace qmcplusplus
{

StressPBC::StressPBC(ParticleSet& ions, ParticleSet& elns, TrialWaveFunction& Psi0, bool firsttime):
  ForceBase(ions, elns), PtclTarg(elns), PtclA(ions), Psi(Psi0), first_time(firsttime)
{
  ReportEngine PRE("StressPBC","StressPBC");
  myName = "StressPBC";
  prefix="StressPBC";
  //Defaults for the chiesa S-wave polynomial filtering.  
  Rcut = 0.4;
  m_exp = 2;
  N_basis = 4;
 // forces = 0.0;
 // forces_ShortRange.resize(Nnuc);
 // forces_ShortRange = 0.0;
 // forces_IonIon=0.0;
  
  //This sets up the long range breakups. 
  
  kcdifferent=false;
  myTableIndex=elns.addTable(PtclA,DT_AOS);
  initBreakup(PtclTarg);
  
  targetconsts=0.0;
  stress_eI_const=0.0;
  stress_ee_const=0.0;
  
  if (firsttime==true)
  { 
	  CalculateIonIonStress();
      firsttime=false;
  }

 // app_log()<<"Ion sr = "<<evaluateSR_AA(PtclA)<< std::endl;
 // app_log()<<"\nIon lr = "<<evaluateLR_AA(PtclA)<< std::endl;
//  app_log()<<"\n Ion const = "<<evalConsts_AA(PtclA)<< std::endl;
 // app_log()<<"\n e-e const = "<<evalConsts_AA(PtclTarg)<< std::endl;
//  app_log()<<"\n e-I const = "<<evalConsts_AB()<< std::endl<< std::endl;
//  evaluateSR_AA();
  
//  stress_IonIon=evaluateSR_AA(PtclA)+evaluateLR_AA(PtclA)+evalConsts_AA(PtclA); //+ evaluateLR_AA(PtclA);
//  stress_eI_const+=evalConsts_AB();
 // stress_ee_const+=evalConsts_AA(PtclTarg);
  
 // app_log()<<"\n====ion-ion stress ====\n"<<stress_IonIon<< std::endl;
 // app_log()<<"\n eI_const = "<<stress_eI_const<< std::endl;
  //app_log()<< "IonIon Force" << std::endl;
 // app_log()<<forces_IonIon<< std::endl; 
 // app_log() << "  Maximum K shell " << AB->MaxKshell << std::endl;
 // app_log() << "  Number of k vectors " << AB->Fk.size() << std::endl;
  
  ///////////////////////////////////////////////////////////////
}

void StressPBC::InitMatrix()
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

void StressPBC::initBreakup(ParticleSet& P)
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
  
  for(int spec=0; spec<NumSpeciesA; spec++)
  {
    
  }
  RealType totQ=0.0;
  for(int iat=0; iat<NptclA; iat++)
    totQ+=Zat[iat] = Zspec[PtclA.GroupID[iat]];
  for(int iat=0; iat<NptclB; iat++)
    totQ+=Qat[iat] = Qspec[P.GroupID[iat]];
 
  kcdifferent = (std::abs(PtclA.Lattice.LR_kc - P.Lattice.LR_kc) > std::numeric_limits<RealType>::epsilon());
  minkc = std::min(PtclA.Lattice.LR_kc,P.Lattice.LR_kc);
  //AB->initBreakup(*PtclB);
  //initBreakup is called only once
  //AB = LRCoulombSingleton::getHandler(*PtclB);
 // AB = LRCoulombSingleton::getDerivHandler(P);
  AA = LRCoulombSingleton::getDerivHandler(P);
  //myConst=evalConsts();
  myRcut=AA->get_rc();//Basis.get_rc();
}






SymTensor<StressPBC::RealType,OHMMS_DIM> StressPBC::evaluateLR_AB(ParticleSet& P)
{
 // const int slab_dir=OHMMS_DIM-1;
  SymTensor<RealType, OHMMS_DIM> res=0.0;
  const StructFact& RhoKA(*(PtclA.SK));
  const StructFact& RhoKB(*(P.SK));

    for(int i=0; i<NumSpeciesA; i++)
    {
      SymTensor<RealType, OHMMS_DIM> esum;
      esum=0.0;
      for(int j=0; j<NumSpeciesB; j++)
      {
#if defined(USE_REAL_STRUCT_FACTOR)
        esum += Qspec[j]*AA->evaluateStress(RhoKA.KLists.kshell
                                      , RhoKA.rhok_r[i],RhoKA.rhok_i[i] , RhoKB.rhok_r[j],RhoKB.rhok_i[j]);
#else
        esum += Qspec[j]*AA->evaluateStress(RhoKA.KLists.kshell, RhoKA.rhok[i],RhoKB.rhok[j]);

#endif
      } //speceln

      res += Zspec[i]*esum;
    }
      

  return res;
  
}

SymTensor<StressPBC::RealType,OHMMS_DIM> StressPBC::evaluateSR_AB(ParticleSet& P)
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
      esum += Qat[jat]*AA->evaluateSR_dstrain(d_ab.dr(nn), d_ab.r(nn));
    }
    //Accumulate pair sums...species charge for atom i.
    res += Zat[iat]*esum;
  }

  return res;
}


SymTensor<StressPBC::RealType,OHMMS_DIM> StressPBC::evaluateSR_AA(ParticleSet& P)
{
  const DistanceTableData &d_aa(*P.DistTables[0]);
  SymTensor<RealType,OHMMS_DIM> stress;
   int NumSpecies = P.getSpeciesSet().TotalNum;
  //RealType res=0.0;
  //Loop over distinct eln-ion pairs
 //app_log()<<"NumSpeciesA = "<<NumSpeciesA<< std::endl;
 
//  int ChargeAttribIndx = P.getSpeciesSet().getAttribute("charge");
 // int MemberAttribIndx =P.getSpeciesSet().getAttribute("membersize");
  
 // std::vector<int> NofSpecies;
 // NofSpecies.resize(NumSpecies);
 // std::vector<int> Zmyspec;
//  Zmyspec.resize(NumSpecies);
//  for(int spec=0; spec<NumSpecies; spec++)
//  {
 //   Zmyspec[spec] = P.getSpeciesSet()(ChargeAttribIndx,spec);
//    NofSpecies[spec] = static_cast<int>(P.getSpeciesSet()(MemberAttribIndx,spec));
//  }
 
 for(int ipart=0; ipart<P.getTotalNum(); ipart++)
  {
    SymTensor<RealType, OHMMS_DIM> esum = 0.0;
    for(int nn=d_aa.M[ipart],jpart=ipart+1; nn<d_aa.M[ipart+1]; nn++,jpart++)
    {
      esum += P.Z[jpart]* AA->evaluateSR_dstrain(d_aa.dr(nn), d_aa.r(nn));
    }
    stress+= P.Z[ipart]*esum;
  }
  
  return stress;
 
}

SymTensor<StressPBC::RealType,OHMMS_DIM> StressPBC::evaluateLR_AA(ParticleSet& P)
{

 // const DistanceTableData &d_aa(*P.DistTables[0]);
  int NumSpecies = P.getSpeciesSet().TotalNum;
  SymTensor<RealType, OHMMS_DIM> stress; 
  const StructFact& PtclRhoK(*(P.SK));
  //vector<TinyVector<RealType,DIM> > grad(P.getTotalNum());
  int ChargeAttribIndx = P.getSpeciesSet().getAttribute("charge");
  int MemberAttribIndx =P.getSpeciesSet().getAttribute("membersize");
  
  std::vector<int> NofSpecies;
  std::vector<int> Zmyspec;
  NofSpecies.resize(NumSpecies);
  Zmyspec.resize(NumSpecies);;

  for(int spec=0; spec<NumSpecies; spec++)
  {
    Zmyspec[spec] = P.getSpeciesSet()(ChargeAttribIndx,spec);
    NofSpecies[spec] = static_cast<int>(P.getSpeciesSet()(MemberAttribIndx,spec));
  }
  

  SymTensor<RealType, OHMMS_DIM> temp;
  for(int spec1=0; spec1<NumSpecies; spec1++)
  {
    RealType Z1 = Zmyspec[spec1];
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
        stress += Z1*Zmyspec[spec2]*temp;
    } //spec2
  }//spec1
  
  return stress;
}


SymTensor<StressPBC::RealType, OHMMS_DIM>
StressPBC::evalConsts_AB()
{
  int nelns = PtclTarg.getTotalNum();
  int nions = PtclA.getTotalNum();

  typedef LRHandlerType::mRealType mRealType;

  SymTensor<mRealType, OHMMS_DIM> Consts=0.0;
  SymTensor<mRealType, OHMMS_DIM> vs_k0 = AA->evaluateSR_k0_dstrain();
  mRealType v1; //single particle energy
  for(int i=0; i<nelns; ++i)
  {
    v1=0.0;
    for(int s=0; s<NumSpeciesA; s++)
      v1 += NofSpeciesA[s]*Zspec[s];
   // v1 *= -.5*Qat[i]*vs_k0;
    //Ve_const(i) = v1;
    Consts += -.5*Qat[i]*vs_k0*v1;
  }
  for(int i=0; i<nions; ++i)
  {
    v1=0.0;
    for(int s=0; s<NumSpeciesB; s++)
      v1 += NofSpeciesB[s]*Qspec[s];
   // v1 *= -.5*Zat[i]*vs_k0;
   // Vi_const(i) = v1;
    Consts += -.5*Zat[i]*vs_k0*v1;
  }

  return SymTensor<RealType,OHMMS_DIM>(Consts);
}

SymTensor<StressPBC::RealType, OHMMS_DIM>
StressPBC::evalConsts_AA(ParticleSet& P)
{
  SymTensor<RealType, OHMMS_DIM> tmpconsts=0.0; // constant term
 
  int NumSpecies = P.getSpeciesSet().TotalNum;
  int NumCenters = P.getTotalNum();
  RealType v1; //single particle energy

  int ChargeAttribIndx = P.getSpeciesSet().getAttribute("charge");
  int MemberAttribIndx =P.getSpeciesSet().getAttribute("membersize");
  
  std::vector<int> NofSpecies;
  std::vector<int> Zmyspec;
  NofSpecies.resize(NumSpecies);
  Zmyspec.resize(NumSpecies);
 
  for(int spec=0; spec<NumSpecies; spec++)
  {
    Zmyspec[spec] = P.getSpeciesSet()(ChargeAttribIndx,spec);
    NofSpecies[spec] = static_cast<int>(P.getSpeciesSet()(MemberAttribIndx,spec));
  }
  
  
//  V_const = 0.0;
  //v_l(r=0) including correction due to the non-periodic direction
  SymTensor<RealType, OHMMS_DIM> vl_r0 = AA->evaluateLR_r0_dstrain();
  //app_log()<<"   PBCAA vl_r0 = "<<vl_r0<< std::endl;
  for(int ipart=0; ipart<NumCenters; ipart++)
  {
  //  v1 =  -.5*Zat[ipart]*Zat[ipart]*vl_r0;
   // V_const(ipart) += v1;
    tmpconsts += -.5*P.Z[ipart]*P.Z[ipart]*vl_r0;
  }
 // if(report)
    app_log() << "   PBCAA self-interaction term \n" << tmpconsts << std::endl;
  //Compute Madelung constant: this is not correct for general cases
  SymTensor<RealType, OHMMS_DIM> MC0 = 0;
  for(int ks=0; ks<AA->Fk.size(); ks++)
    MC0 += AA->dFk_dstrain[ks];
  MC0 = 0.5*(MC0 - vl_r0);
  //Neutraling background term
  SymTensor<RealType, OHMMS_DIM> vs_k0=AA->evaluateSR_k0_dstrain(); //v_s(k=0)
  for(int ipart=0; ipart<NumCenters; ipart++)
  {
    v1 = 0.0;
    for(int spec=0; spec<NumSpecies; spec++)
      v1 += -.5*P.Z[ipart]*NofSpecies[spec]*Zmyspec[spec];
   // v1 *= -.5*Zat[ipart]*vs_k0;
 //   V_const(ipart) += v1;
    tmpconsts += v1*vs_k0;
  }
 // if(report)
    app_log() << "   PBCAA total constant \n" << tmpconsts << std::endl;
  //app_log() << "   MC0 of PBCAA " << MC0 << std::endl;
  return tmpconsts;
}



StressPBC::Return_t
StressPBC::evaluate(ParticleSet& P)
{
  stress=0.0;
  
  stress+=evaluateLR_AB(PtclTarg);
  stress+=evaluateSR_AB(PtclTarg);
  stress+=stress_eI_const;
  
  stress+=evaluateLR_AA(PtclTarg);
  stress+= evaluateSR_AA(PtclTarg);
  stress+=stress_ee_const;
  
  stress+=stress_IonIon;

  stress+=evaluateKineticSymTensor(P);
 
  //stress/=(-1.0*P.Lattice.Volume);
  const RealType vinv(-1.0/P.Lattice.Volume);
  stress *= vinv;

  return 0.0;
}

SymTensor<StressPBC::RealType,OHMMS_DIM> StressPBC::evaluateKineticSymTensor(ParticleSet& P)
{
  OrbitalBase::HessVector_t grad_grad_psi;
  Psi.evaluateHessian(P,grad_grad_psi);
  SymTensor<RealType,OHMMS_DIM> kinetic_tensor;
  Tensor<ComplexType, OHMMS_DIM> complex_ktensor;

  
  for (int iat=0; iat<P.getTotalNum(); iat++)
  {
    const RealType minv(1.0/P.Mass[iat]);
    complex_ktensor+=outerProduct(P.G[iat],P.G[iat])*static_cast<ParticleSet::ParticleValue_t>(minv);
    complex_ktensor+=grad_grad_psi[iat]*minv;
    //complex_ktensor+=(grad_grad_psi[iat] + outerProduct(P.G[iat],P.G[iat]))*(1.0/(P.Mass[iat]));
  }
  
  for (int i=0; i<OHMMS_DIM; i++)
	for(int j=i; j<OHMMS_DIM; j++)
	{
	  kinetic_tensor(i,j)=complex_ktensor(i,j).real();
	}
  return kinetic_tensor;

}

StressPBC::Return_t StressPBC::g_filter(RealType r)
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

bool StressPBC::put(xmlNodePtr cur)
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
  app_log() <<"    StressPBC Parameters"<< std::endl;
  app_log() <<"        StressPBC::Rcut="<<Rcut<< std::endl;
  app_log() <<"        StressPBC::N_basis="<<N_basis<< std::endl;
  app_log() <<"        StressPBC::m_exp="<<m_exp<< std::endl;
  InitMatrix();
  return true;
}

QMCHamiltonianBase* StressPBC::makeClone(ParticleSet& qp, TrialWaveFunction& psi)
{
  StressPBC* tmp = new StressPBC(PtclA, qp, psi, false);
  tmp->targetconsts=targetconsts;
  tmp->stress_IonIon=stress_IonIon;
  tmp->stress_ee_const=stress_ee_const;
  tmp->stress_eI_const=stress_eI_const;
  
  return tmp;
}
}


