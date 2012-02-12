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
#include "qmc_common.h"
#include <QMCApp/ParticleSetPool.h>

namespace qmcplusplus 
{
  ///initialize static data
  ParticleSetPool* qmc_common::ptcl_pool=0;

  ParticleSetPool* qmc_common::getParticleSetPool(Communicate* mycomm)
  {
    if(ptcl_pool==0)
      ptcl_pool=new ParticleSetPool(mycomm);
    return ptcl_pool;
  }
}
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 5388 $   $Date: 2011-12-02 08:45:44 -0500 (Fri, 02 Dec 2011) $
 * $Id: Configuration.h 5388 2011-12-02 13:45:44Z jnkim $ 
 ***************************************************************************/
