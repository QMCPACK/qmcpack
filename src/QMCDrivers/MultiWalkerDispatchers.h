//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_MWDISPATCHERS_H
#define QMCPLUSPLUS_MWDISPATCHERS_H

#include <TWFdispatcher.h>
#include <Hdispatcher.h>

namespace qmcplusplus
{
// forward declaration
class TWFdispatcher;
class Hdispatcher;

class MultiWalkerDispatchers
{
public:

  MultiWalkerDispatchers(bool use_batch) : twf_dispatcher_(use_batch), ham_dispatcher_(use_batch) {}

  const TWFdispatcher twf_dispatcher_;
  const Hdispatcher ham_dispatcher_;
};
}
#endif
