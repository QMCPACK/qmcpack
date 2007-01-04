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
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#include "Particle/HDFWalkerInputManager.h"
#include "Message/Communicate.h"
#include "OhmmsData/AttributeSet.h"
#include "Particle/HDFWalkerInput0.h"
#include "Particle/HDFWalkerInputCollect.h"

namespace qmcplusplus {

  HDFWalkerInputManager::HDFWalkerInputManager(MCWalkerConfiguration& w):
  Wref(w), CurrentFileRoot("invalid") {
  }
  HDFWalkerInputManager::~HDFWalkerInputManager() {
  }

  bool HDFWalkerInputManager::put(xmlNodePtr cur) {

    app_error() << "HDFWalkerInputManager::put(xmlNodePtr cur) "
      << " is not implemented." << endl;
    return false;
  }

  bool HDFWalkerInputManager::put(std::vector<xmlNodePtr> wset) {
    int pid=OHMMS::Controller->mycontext(); 
    int nfile=wset.size();
    for(int ifile=0; ifile<nfile; ifile++) {
      string cfile("invalid"), target("e"), collect("no");
      int anode=-1, nwalkers=-1, nblocks=1;
      OhmmsAttributeSet pAttrib;
      pAttrib.add(cfile,"href"); pAttrib.add(cfile,"file"); 
      pAttrib.add(target,"target"); pAttrib.add(target,"ref"); 
      pAttrib.add(anode,"node");
      pAttrib.add(nwalkers,"walkers");
      pAttrib.add(nblocks,"rewind");
      pAttrib.add(collect,"collect");
      pAttrib.put(wset[ifile]);

      if(collect == "no") { // use old method
        int pid_target= (anode<0)?pid:anode;
        if(pid_target == pid && cfile != "invalid") {
          HDFWalkerInput0 WO(cfile); 
          WO.append(Wref,nblocks);
          //WO.append(Wref,nwalkers);
          //read random state
          WO.getRandomState(true);
          CurrentFileRoot = cfile;
        }
      } else {
        HDFWalkerInputCollect WO(cfile);
        WO.put(Wref,nblocks);
        CurrentFileRoot = cfile;
      }
    }
    return true;
  }

  void HDFWalkerInputManager::rewind(const std::string& h5root, int blocks) {
    HDFWalkerInputCollect WO(h5root);
    WO.rewind(Wref,blocks);
  }
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
