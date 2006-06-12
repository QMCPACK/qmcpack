//////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Jeongnim Kim
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
/** @file ECPComponentBuilderBuilder.h
 * @brief Declaration of a builder class for an ECP component for an ionic type
 */
#ifndef QMCPLUSPLUS_ECPCOMPONENT_BUILDER_H
#define QMCPLUSPLUS_ECPCOMPONENT_BUILDER_H
#include "Particle/DistanceTableData.h"
#include "QMCHamiltonians/LocalECPotential.h"
#include "QMCHamiltonians/NonLocalECPotential.h"

namespace qmcplusplus {

  struct ECPComponentBuilder: public QMCTraits {

    typedef LocalECPotential::GridType GridType;
    typedef LocalECPotential::RadialPotentialType RadialPotentialType;

    int NumNonLocal;
    int Lmax;
    RealType Zeff;
    string Species;
    RadialPotentialType* pp_loc;
    NonLocalECPComponent* pp_nonloc;
    map<string,int> angMon;

    ECPComponentBuilder(const string& aname);

    bool parse(const string& fname);
    bool put(xmlNodePtr cur);
    void addSemiLocal(xmlNodePtr cur);
    void buildLocal(xmlNodePtr cur);
    void buildSemiLocalAndLocal(vector<xmlNodePtr>& semiPtr);
    RadialPotentialType* createVr(xmlNodePtr cur, GridType* agrid);
  };

}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
