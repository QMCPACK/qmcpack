//////////////////////////////////////////////////////////////////
// (c) Copyright 2005- by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim
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
#include "QMCDrivers/WalkerControlBase.h"
namespace qmcplusplus {

  int WalkerControlBase::branch(int iter, MCWalkerConfiguration& W, RealType trigger) {

    sortWalkers(W);

    int nw_tot = copyWalkers(W);

    //set Weight and Multiplicity to default values
    MCWalkerConfiguration::iterator it(W.begin()),it_end(W.end());
    while(it != it_end) {
      (*it)->Weight= 1.0;
      (*it)->Multiplicity=1.0;
      ++it;
    }

    return nw_tot;
  }

  void WalkerControlBase::sortWalkers(MCWalkerConfiguration& W) {

    MCWalkerConfiguration::iterator it(W.begin());

    vector<Walker_t*> bad;
    NumWalkers=0;
    MCWalkerConfiguration::iterator it_end(W.end());
    while(it != it_end) {
      int nc = std::min(static_cast<int>((*it)->Multiplicity),MaxCopy);
      if(nc) {
        NumWalkers += nc;
        good_w.push_back(*it);
        ncopy_w.push_back(nc-1);
      } else {
        bad.push_back(*it);
      }
      ++it;
    }

    //remove bad walkers empty the container
    for(int i=0; i<bad.size(); i++) delete bad[i];

    if(good_w.empty()) {
      app_error() << "All the walkers have died. Abort. " << endl;
      OHMMS::Controller->abort();
    }

    int sizeofgood = good_w.size();

    //check if the projected number of walkers is too small or too large
    if(NumWalkers>Nmax) {
      int nsub=0;
      int nsub_target=NumWalkers-static_cast<int>(0.9*Nmax);
      int i=0;
      while(i< sizeofgood && nsub<nsub_target) {
        if(ncopy_w[i]) {ncopy_w[i]--; nsub++;}
        ++i;
      }
      NumWalkers -= nsub;
    } else  if(NumWalkers < Nmin) {
      int nadd=0;
      int nadd_target = static_cast<int>(Nmin*1.1)-NumWalkers;
      if(nadd_target> sizeofgood) {
        app_warning() << "The number of walkers is running low. Requested walkers " 
          << nadd_target << " good walkers = " << sizeofgood << endl;
      }
      int i=0;
      while(i< sizeofgood && nadd<nadd_target) {
        ncopy_w[i]++; ++nadd;++i;
      }
      NumWalkers +=  nadd;
    }
  }

  int WalkerControlBase::copyWalkers(MCWalkerConfiguration& W) {
    //clear the WalkerList to populate them with the good walkers
    W.clear();
    W.insert(W.begin(), good_w.begin(), good_w.end());

    int cur_walker = good_w.size();
    for(int i=0; i<good_w.size(); i++) { //,ie+=ncols) {
      for(int j=0; j<ncopy_w[i]; j++, cur_walker++) {
        W.push_back(new Walker_t(*(good_w[i])));
      }
    }

    //clear good_w and ncopy_w for the next branch
    good_w.clear();
    ncopy_w.clear();
    return W.getActiveWalkers();
  }

//#if defined(HAVE_MPI)
//  WalkerControlBase::RealType WalkerControlBase::average(RealType eavg, RealType wgt) {
//    gEavgWgt[0]=eavg;
//    gEavgWgt[1]=wgt;
//    gsum(gEavgWgt,0);
//    return gEavgWgt[0]/gEavgWgt[1];
//  }
//#endif

}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

