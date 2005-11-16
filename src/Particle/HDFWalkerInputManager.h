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
#ifndef QMCPLUSPLUS_WALKER_INPUT_MANAGER_H
#define QMCPLUSPLUS_WALKER_INPUT_MANAGER_H

#include "Particle/MCWalkerConfiguration.h"
#include "OhmmsData/HDFAttribIO.h"
#include <queue>

namespace qmcplusplus {

  class HDFWalkerInputManager {

    MCWalkerConfiguration& Wref;
    std::string CurrentFileRoot;

    public:

      HDFWalkerInputManager(MCWalkerConfiguration& w);
      ~HDFWalkerInputManager();
      bool put(xmlNodePtr cur);
      bool put(std::vector<xmlNodePtr> mset);
      std::string getLastFile() { return CurrentFileRoot;}

      void rewind(const std::string& h5root, int blocks);
    
  };
}

#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
