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
#include "OhmmsData/ParameterSet.h"
#include "QMCDrivers/DMC/WalkerControlFactory.h"
#include "QMCDrivers/DMC/WalkerReconfiguration.h"
#if defined(HAVE_MPI)
#include "QMCDrivers/DMC/WalkerControlMPI.h"
#include "QMCDrivers/DMC/WalkerReconfigurationMPI.h"
#endif
namespace qmcplusplus {

#if defined(HAVE_MPI)
  WalkerControlBase* CreateWalkerController(
      bool reconfig, int& swapmode, int nideal,
      int nmax, int nmin, WalkerControlBase* wc,
      Communicate* comm) {

      int ncontexts = comm->ncontexts();

      //overwrite the SwapMode for a single-node run
      if(ncontexts == 1) {swapmode=0;}

      if(reconfig) {
        nmax=nideal/ncontexts;
        nmin=nideal/ncontexts;
      } else {
        if(swapmode) {
          int npernode=nideal/ncontexts;
          nmax=2*npernode+1;
          nmin=npernode/5+1;
        } else {
          nmax=2*nideal;
          nmin=nideal/2;
        }
      }

      if(wc) {
        if(swapmode != wc->SwapMode) {
          delete wc;
          wc=0;
        }
      } 

      if(wc == 0) {
        if(swapmode) {
          if(reconfig)  {
            app_log() << "  Using WalkerReconfigurationMPI for population control." << endl;
            wc = new WalkerReconfigurationMPI(comm);
          } else {
            app_log() << "  Using WalkerControlMPI for dynamic population control." << endl;
            wc = new WalkerControlMPI(comm);
          }
        } else {
          if(reconfig)  {
            app_log() << "  Using WalkerReconfiguration for population control." << endl;
            wc = new WalkerReconfiguration;
          } else {
            app_log() << "  Using WalkerControlBase for dynamic population control." << endl;
            wc = new WalkerControlBase;
          }
        }
      }

      wc->Nmin=nmin;
      wc->Nmax=nmax;
      return wc;
    }
#else
  WalkerControlBase* CreateWalkerController(
      bool reconfig, int& swapmode, int nideal,
      int nmax, int nmin, WalkerControlBase* wc,
      Communicate* comm) {
    //reset to 0 so that never ask the same question
    swapmode = 0;
    //if(nmax<0) nmax=2*nideal;
    //if(nmin<0) nmin=nideal/2;

    //if(wc== 0) wc= new WalkerControlBase;
    if(wc== 0) {
      if(reconfig) {
        app_log() << "  Using a fixed number of walkers by reconfiguration." << endl;
        wc = new WalkerReconfiguration(comm);
        wc->Nmax=nideal;
        wc->Nmin=nideal;
      } else {
        app_log() << "  Using a WalkerControlBase with population fluctations." << endl;
        wc = new WalkerControlBase(comm);
        wc->Nmax=2*nideal;
        wc->Nmin=nideal/2;
      }
    }
    return wc;
  }

  WalkerControlBase* createWalkerController(int nwtot, Communicate* comm, xmlNodePtr cur) 
  {

    app_log() << "  Creating WalkerController: current number of walkers = " << nwtot << endl;

    ///set of parameters
    int nmax=0;
    int nmin=0;
    string reconfig("no");
    ParameterSet m_param;
    m_param.add(nwtot,"targetWalkers","int"); 
    m_param.add(nwtot,"targetwalkers","int"); 
    m_param.add(reconfig,"reconfiguration","string");
    m_param.put(cur);

    //if(nmax<0) nmax=2*nideal;
    //if(nmin<0) nmin=nideal/2;

    WalkerControlBase* wc=0;

    //if(wc== 0) wc= new WalkerControlBase;
    if(reconfig == "yes") {
      app_log() << "  Using a fixed number of walkers by reconfiguration." << endl;
      wc = new WalkerReconfiguration(comm);
      wc->Nmax=nwtot;
      wc->Nmin=nwtot;
    } else {
      app_log() << "  Using a WalkerControlBase with population fluctations." << endl;
      app_log() << "  Target number of walkers = " << nwtot << endl;
      wc = new WalkerControlBase(comm);
      wc->Nmax=2*nwtot;
      wc->Nmin=nwtot/2;
    }

    return wc;
  }
#endif
}
/***************************************************************************
 * $RCSfile: WalkerControlFactory.cpp,v $   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

