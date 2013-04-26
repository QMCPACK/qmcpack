//////////////////////////////////////////////////////////////////
// (c) Copyright 2006-  by Kris Delaney and Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#include <LongRange/EwaldHandler.h>

namespace qmcplusplus
{
void EwaldHandler::initBreakup(ParticleSet& ref)
{
  SuperCellEnum=ref.Lattice.SuperCellEnum;
  LR_rc=ref.Lattice.LR_rc;
  LR_kc=ref.Lattice.LR_kc;
  Sigma=3.5;
  //determine the sigma
  while(erfc(Sigma)/LR_rc>1e-10)
  {
    Sigma+=0.1;
  }
  app_log() << "   EwaldHandler Sigma/LR_rc = " << Sigma ;
  Sigma/=ref.Lattice.LR_rc;
  app_log() << "  Sigma=" << Sigma  << endl;
  Volume=ref.Lattice.Volume;
  PreFactors=0.0;
  //See A.18
  if(SuperCellEnum == SUPERCELL_SLAB)
  {
    Area=abs(ref.Lattice.R(0,0)*ref.Lattice.R(1,1)-ref.Lattice.R(0,1)*ref.Lattice.R(1,0));
    PreFactors[0]=2.0*M_PI/Area; //\f$ \frac{2\pi}{A}\f$
    PreFactors[1]=2.0*std::sqrt(M_PI)/(Sigma*Area);//\f$  \frac{2\pi}{A}\frac{1}{Sigma\pi}\f$
    PreFactors[2]=1.0/(2*Sigma);//Used for the k-dependent term
  }
  fillFk(ref.SK->KLists);
}

EwaldHandler::EwaldHandler(const EwaldHandler& aLR, ParticleSet& ref):
  LRHandlerBase(aLR), Sigma(aLR.Sigma), Volume(aLR.Volume), Area(aLR.Area)
  , PreFactors(aLR.PreFactors)
{
  SuperCellEnum = aLR.SuperCellEnum;
  if(SuperCellEnum== SUPERCELL_SLAB)
    kMag=aLR.kMag;
}

void EwaldHandler::fillFk(KContainer& KList)
{
  Fk.resize(KList.kpts_cart.size());
  const vector<int>& kshell(KList.kshell);
  MaxKshell=kshell.size()-1;
  Fk_symm.resize(MaxKshell);
  kMag.resize(MaxKshell);
  if(SuperCellEnum==SUPERCELL_SLAB)
  {
    RealType knorm=M_PI/Area;
    RealType ksum=0.0;
    RealType oneovertwosigma=1.0/(2.0*Sigma);
    for(int ks=0,ki=0; ks<Fk_symm.size(); ks++)
    {
      kMag[ks]=std::sqrt(KList.ksq[ki]);
      RealType uk=knorm/kMag[ks];//pi/(A*k)
      Fk_symm[ks]=uk;
      while(ki<KList.kshell[ks+1] && ki<Fk.size())
        Fk[ki++]=uk;
      ksum += uk*erfc(kMag[ks]*oneovertwosigma)*(KList.kshell[ks+1]-KList.kshell[ks]);
    }
    //correction for the self-energy returned by evaluateLR_r0
    PreFactors[3]=2.0*(std::sqrt(M_PI)/(Area*Sigma)-ksum);
  }
  else
  {
#if OHMMS_DIM==2
    RealType kgauss=1.0/(4*Sigma*Sigma);
    RealType knorm=2*M_PI/Volume;
    const RealType acclog=std::abs(std::log(1.0e-10));
    for(int ks=0,ki=0; ks<Fk_symm.size(); ks++)
    {
      RealType t2e=KList.ksq[ki]*kgauss;
      RealType uk=knorm*std::exp(-t2e)/KList.ksq[ki];
      Fk_symm[ks]=uk;
      while(ki<KList.kshell[ks+1] && ki<Fk.size())
        Fk[ki++]=uk;
    }
    PreFactors[3]=0.0;
#elif OHMMS_DIM==3
    RealType kgauss=1.0/(4*Sigma*Sigma);
    RealType knorm=4*M_PI/Volume;
    const RealType acclog=std::abs(std::log(1.0e-10));
    for(int ks=0,ki=0; ks<Fk_symm.size(); ks++)
    {
      RealType t2e=KList.ksq[ki]*kgauss;
      RealType uk=knorm*std::exp(-t2e)/KList.ksq[ki];
      Fk_symm[ks]=uk;
      while(ki<KList.kshell[ks+1] && ki<Fk.size())
        Fk[ki++]=uk;
    }
    PreFactors[3]=0.0;
#endif
  }
  app_log().flush();
}

EwaldHandler::RealType
EwaldHandler::evaluate_slab(RealType z, const vector<int>& kshell
                            , const ComplexType* restrict eikr_i, const ComplexType* restrict eikr_j)
{
  RealType zp=z*Sigma;
  RealType vk=-SlabFunc0(z,zp);
  //cout << "### SLAB " << z << " " << zp << endl;
  for(int ks=0,ki=0; ks<MaxKshell; ks++)
  {
    RealType u=0;//\sum Real (e^ikr_i e^(-ikr_j))
    for(; ki<kshell[ks+1]; ki++,eikr_i++,eikr_j++)
      u += ((*eikr_i).real()*(*eikr_j).real()+(*eikr_i).imag()*(*eikr_j).imag());
    vk += u*Fk_symm[ks]*SlabFuncK(ks,z,zp);
  }
  return vk;
}
}
/***************************************************************************
 * $RCSfile$   $Author: bkclark $
 * $Revision: 3798 $   $Date: 2009-04-29 00:38:29 -0500 (Wed, 29 Apr 2009) $
 * $Id: EwaldHandler.cpp 3798 2009-04-29 05:38:29Z bkclark $
 ***************************************************************************/
