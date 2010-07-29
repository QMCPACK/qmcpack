//////////////////////////////////////////////////////////////////
// (c) Copyright 2003 by Jeongnim Kim and Jordan Vincent
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
#include "QMCWaveFunctions/MolecularOrbitals/GTOBuilder.h"
using namespace std;
namespace qmcplusplus {

  GTOBuilder::GTOBuilder(xmlNodePtr cur): Normalized(true), m_orbitals(0) {
    if(cur != NULL) {
      putCommon(cur);
    }
  }

  bool GTOBuilder::putCommon(xmlNodePtr cur) {
    const xmlChar* a=xmlGetProp(cur,(const xmlChar*)"normalized");
    if(a) {
      if(xmlStrEqual(a,(const xmlChar*)"no")) Normalized=false;
    }
    return true;
  }

  bool
  GTOBuilder::addRadialOrbital(xmlNodePtr cur, const QuantumNumberType& nlms) { 

    if(!m_orbitals) {
      ERRORMSG("m_orbitals, SphericalOrbitals<ROT,GT>*, is not initialized")
      return false;
    }

    RadialOrbitalType* radorb= new RadialOrbitalType(nlms[q_l],Normalized);
    radorb->putBasisGroup(cur);
    m_orbitals->Rnl.push_back(radorb);
    m_orbitals->RnlID.push_back(nlms);

    return true;
  }

}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
