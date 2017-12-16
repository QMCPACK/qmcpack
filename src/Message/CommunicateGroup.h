//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
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
