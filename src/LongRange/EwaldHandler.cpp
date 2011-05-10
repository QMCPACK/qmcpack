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
    LR_rc=ref.Lattice.LR_rc;
    LR_kc=ref.Lattice.LR_kc;
    Sigma=3.5;
    while(erfc(Sigma)/LR_rc>1e-10)
    {
      Sigma+=0.1;
    }
    app_log() << "   EwaldHandler Sigma/LR_rc = " << Sigma << endl;

    Sigma/=ref.Lattice.LR_rc;
    Volume=ref.Lattice.Volume;
    fillFk(ref.SK->KLists);
  }

  EwaldHandler::EwaldHandler(const EwaldHandler& aLR, ParticleSet& ref):
    LRHandlerBase(aLR), Sigma(aLR.Sigma), Volume(aLR.Volume)
  {
    //nothing to be done
  }

  void EwaldHandler::fillFk(KContainer& KList)
  {
    RealType kgauss=1.0/(4*Sigma*Sigma);
    RealType knorm=4*M_PI/Volume;
    const RealType acclog=std::abs(std::log(1.0e-10));

    Fk.resize(KList.kpts_cart.size());
    const vector<int>& kshell(KList.kshell);
    if(MaxKshell >= kshell.size()) MaxKshell=kshell.size()-1;
    Fk_symm.resize(MaxKshell);
    int ks_max=0;
    for(int ks=0,ki=0; ks<Fk_symm.size(); ks++)
    {

      RealType t2e=KList.ksq[ki]*kgauss;
      if(-t2e<acclog) ks_max=ks;
      RealType uk=knorm*std::exp(-t2e)/KList.ksq[ki];
      Fk_symm[ks]=uk;
      while(ki<KList.kshell[ks+1] && ki<Fk.size()) Fk[ki++]=uk;
    }
    app_log() <<  "  EwaldHanlder MaxKshell = " << ks_max << endl;

    app_log().flush();
  }
}
/***************************************************************************
 * $RCSfile$   $Author: bkclark $
 * $Revision: 3798 $   $Date: 2009-04-29 00:38:29 -0500 (Wed, 29 Apr 2009) $
 * $Id: EwaldHandler.cpp 3798 2009-04-29 05:38:29Z bkclark $
 ***************************************************************************/
