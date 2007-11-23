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
#ifndef QMCPLUSPLUS_WALKER_INPUT_MANAGER_H
#define QMCPLUSPLUS_WALKER_INPUT_MANAGER_H

#include "Particle/MCWalkerConfiguration.h"
#include "OhmmsData/HDFAttribIO.h"
#include <queue>

class Communicate;

namespace qmcplusplus {

  class HDFWalkerInputManager {

    MCWalkerConfiguration& targetW;
    Communicate* myComm;
    std::string CurrentFileRoot;

    public:

      HDFWalkerInputManager(MCWalkerConfiguration& w, Communicate* c);
      ~HDFWalkerInputManager();
      bool put(xmlNodePtr cur);
      //bool put(std::vector<xmlNodePtr>& mset, int pid);
      //bool put(std::vector<xmlNodePtr>& mset, Communicate* comm);
      std::string getFileRoot() { return CurrentFileRoot;}

      void rewind(const std::string& h5root, int blocks);
    
  };
}

#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
