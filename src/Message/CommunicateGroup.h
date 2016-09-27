//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Inc.
//                    Jeremy McMinnis, jmcminis@gmail.com, Navar Inc.
//
// File created by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Inc.
//////////////////////////////////////////////////////////////////////////////////////
    
    



#ifndef QMCPACK_COMMUNICATIONGROUP_H
#define QMCPACK_COMMUNICATIONGROUP_H

#include "Message/Communicate.h"

struct CommunicateGroup: public Communicate
{
  ///parent communicator
  Communicate& parent;
  ///constructor
  CommunicateGroup(Communicate& acomm, int ndiv=1);
};
#endif //QMCPACK_COMMUNICATIONGROUP_H 
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/
