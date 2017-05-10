//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



#include <LongRange/TwoDEwaldHandler.h>

namespace qmcplusplus
{
void TwoDEwaldHandler::initBreakup(ParticleSet& ref)
{
  LR_rc=ref.Lattice.LR_rc;
  LR_kc=ref.Lattice.LR_kc;
  Sigma=LR_kc;
  //determine the sigma
  while(erfc(Sigma)/LR_rc>1e-16)
  {
    Sigma+=0.1;
  }
  app_log() << "   TwoDEwaldHandler Sigma/LR_rc = " << Sigma ;
  Sigma/=LR_rc;
  app_log() << "  Sigma=" << Sigma  << std::endl;
  Volume=ref.Lattice.Volume;
  fillFk(ref.SK->KLists);
}

TwoDEwaldHandler::TwoDEwaldHandler(const TwoDEwaldHandler& aLR, ParticleSet& ref):
  LRHandlerBase(aLR), Sigma(aLR.Sigma), Volume(aLR.Volume)
{
}

void TwoDEwaldHandler::fillFk(KContainer& KList)
{
  Fk.resize(KList.kpts_cart.size());
  const std::vector<int>& kshell(KList.kshell);
  MaxKshell=kshell.size()-1;
  Fk_symm.resize(MaxKshell);
  kMag.resize(MaxKshell);
  RealType knorm=2.0*M_PI/Volume;
  RealType oneovertwosigma=1.0/(2.0*Sigma);
  for(int ks=0,ki=0; ks<Fk_symm.size(); ks++)
  {
    kMag[ks]=std::sqrt(KList.ksq[ki]);
    RealType uk=knorm*erfc(kMag[ks]*oneovertwosigma)/kMag[ks];
    Fk_symm[ks]=uk;
    while(ki<KList.kshell[ks+1] && ki<Fk.size())
      Fk[ki++]=uk;
//       app_log()<<kMag[ks]<<" "<<uk<< std::endl;
  }
  app_log().flush();
}
}
