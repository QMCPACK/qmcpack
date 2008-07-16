//////////////////////////////////////////////////////////////////
// (c) Copyright 2006-  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
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
/** @file LRHandlerTemp.h
 * @brief Define a LRHandler with two template parameters
 */
#include "LongRange/LRCoulombSingleton.h"
#include <numeric>

namespace qmcplusplus {
  //initialization of the static data
  LRCoulombSingleton::LRHandlerType* LRCoulombSingleton::CoulombHandler=0;

  LRCoulombSingleton::LRHandlerType*
    LRCoulombSingleton::getHandler(ParticleSet& ref) {
      if(CoulombHandler ==0) {
        app_log() << "  Create CoulombHandler. " << endl;
        CoulombHandler=new LRHandlerType(ref);
        CoulombHandler->initBreakup(ref);
        return CoulombHandler;
      }
      else
      {
        app_log() << "  Copy CoulombHandler. " << endl;
        return new LRHandlerType(*CoulombHandler,ref);
      }
    }

  LRCoulombSingleton::RadFunctorType*
    LRCoulombSingleton::createSpline4RbyVs(LRHandlerType* aLR, RealType rcut,
        GridType* agrid)
    {
      if(agrid == 0) 
      {
        agrid = new GridType;
        agrid->set(0.0,rcut,1001);
      }

      int ng=agrid->size();
      vector<RealType> v(ng);
      RealType r=(*agrid)[0];

      //check if the first point is not zero
      v[0]=(r>numeric_limits<RealType>::epsilon())? r*aLR->evaluate(r,1.0/r):0.0; 
      for(int ig=1; ig<ng-1; ig++) {
        r=(*agrid)[ig];
        v[ig]=r*aLR->evaluate(r,1.0/r);
      }
      v[0] = 2.0*v[1] - v[2];

      v[ng-1]=0.0;
      RadFunctorType* V0=new RadFunctorType(agrid,v);
      RealType deriv=(v[1]-v[0])/((*agrid)[1]-(*agrid)[0]);
      V0->spline(0,deriv,ng-1,0.0);

      return V0;
    }
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/
