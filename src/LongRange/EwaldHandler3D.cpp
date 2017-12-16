//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//
// File created by: Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
    
    



#include <LongRange/EwaldHandler3D.h>

namespace qmcplusplus
{
void EwaldHandler3D::initBreakup(ParticleSet& ref)
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
  app_log() << "   EwaldHandler3D Sigma/LR_rc = " << Sigma ;
  Sigma/=ref.Lattice.LR_rc;
  app_log() << "  Sigma=" << Sigma  << endl;
  Volume=ref.Lattice.Volume;
  PreFactors=0.0;
  //See A.18
//  if(SuperCellEnum == SUPERCELL_SLAB)
//  {
//    Area=abs(ref.Lattice.R(0,0)*ref.Lattice.R(1,1)-ref.Lattice.R(0,1)*ref.Lattice.R(1,0));
//    PreFactors[0]=2.0*M_PI/Area; //\f$ \frac{2\pi}{A}\f$
//    PreFactors[1]=2.0*std::sqrt(M_PI)/(Sigma*Area);//\f$  \frac{2\pi}{A}\frac{1}{Sigma\pi}\f$
//    PreFactors[2]=1.0/(2*Sigma);//Used for the k-dependent term
//  }
  app_log()<< "Made it to fill Fk\n";
  fillFk(ref.SK->KLists);
  fillYkgstrain(ref.SK->KLists);
  filldFk_dk(ref.SK->KLists);
  app_log()<< "Done with fk\n";
}

EwaldHandler3D::EwaldHandler3D(const EwaldHandler3D& aLR, ParticleSet& ref):
  LRHandlerBase(aLR), Sigma(aLR.Sigma), Volume(aLR.Volume), Area(aLR.Area)
  , PreFactors(aLR.PreFactors)
{
  SuperCellEnum = aLR.SuperCellEnum;
//  if(SuperCellEnum== SUPERCELL_SLAB)
//    kMag=aLR.kMag;
}

void EwaldHandler3D::fillFk(KContainer& KList)
{
  Fk.resize(KList.kpts_cart.size());
  const vector<int>& kshell(KList.kshell);
  MaxKshell=kshell.size()-1;
  Fk_symm.resize(MaxKshell);
  kMag.resize(MaxKshell);
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
  app_log().flush();
}
}
