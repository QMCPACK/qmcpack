//////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Jeongnim Kim
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
/** @file qmc_common.h
 *
 * Declare static data that persists throughout simulation
 */
#ifndef QMCPLUSPLUS_GLOBAL_OBJECTS_H
#define QMCPLUSPLUS_GLOBAL_OBJECTS_H

#include <Configuration.h>

namespace qmcplusplus 
{
  class ParticleSetPool;

  struct qmc_common
  {
    ///@ global objects
    static ParticleSetPool* ptcl_pool;
    ///@}

    ///return ptcl_pool 
    static ParticleSetPool* getParticleSetPool(Communicate* mycomm);
  };
}

#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 5388 $   $Date: 2011-12-02 08:45:44 -0500 (Fri, 02 Dec 2011) $
 * $Id: Configuration.h 5388 2011-12-02 13:45:44Z jnkim $ 
 ***************************************************************************/
