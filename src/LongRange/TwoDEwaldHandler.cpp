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
  app_log() << "  Sigma=" << Sigma  << endl;
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
  const vector<int>& kshell(KList.kshell);
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
//       app_log()<<kMag[ks]<<" "<<uk<<endl;
  }
  app_log().flush();
}
}
/***************************************************************************
 * $RCSfile$   $Author: bkclark $
 * $Revision: 3798 $   $Date: 2009-04-29 00:38:29 -0500 (Wed, 29 Apr 2009) $
 * $Id: TwoDEwaldHandler.cpp 3798 2009-04-29 05:38:29Z bkclark $
 ***************************************************************************/
