//////////////////////////////////////////////////////////////////
// (c) Copyright 2007-  by Ken Esler and Jeongnim Kim
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
#ifndef QMCPLUSPLUS_EINSPLINE_SET_BUILDER_H
#define QMCPLUSPLUS_EINSPLINE_SET_BUILDER_H

#include "QMCWaveFunctions/BasisSetBase.h"
#include "QMCWaveFunctions/EinsplineSet.h"
#include "Numerics/HDFNumericAttrib.h"

namespace qmcplusplus {
  
  class EinsplineSetBuilder : public BasisSetBuilder {
  public:
    //////////////////////
    // Type definitions //
    //////////////////////
    typedef map<string,ParticleSet*> PtclPoolType;

    EinsplineSetBuilder(ParticleSet& p, PtclPoolType& psets, xmlNodePtr cur);
    
    ~EinsplineSetBuilder();
    
    bool put (xmlNodePtr cur);

        /** initialize the Antisymmetric wave function for electrons
     *@param cur the current xml node
     */
    SPOSetBase* createSPOSet(xmlNodePtr cur);
    
  protected:
    

  };

}


#endif
