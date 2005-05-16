//////////////////////////////////////////////////////////////////
// (c) Copyright 2003- by Jeongnim Kim
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
//   Department of Physics, Ohio State University
//   Ohio Supercomputer Center
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
/**@file WaveFunctionPool.cpp
 * @brief Implements WaveFunctionPool operators.
 */
#include "QMCApp/WaveFunctionPool.h"
#include "QMCApp/ParticleSetPool.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "QMCWaveFunctions/JastrowBuilder.h"
#include "QMCWaveFunctions/MolecularOrbitals/MolecularOrbitalBuilder.h"
using namespace std;
#include "OhmmsData/AttributeSet.h"
#include "Utilities/OhmmsInfo.h"

namespace ohmmsqmc {
  
  WaveFunctionPool::WaveFunctionPool(const char* aname):OhmmsElementBase(aname){ }

  bool WaveFunctionPool::put(xmlNodePtr cur) {

    string id("psi0"), target("e"), role("extra");
    OhmmsAttributeSet pAttrib;
    pAttrib.add(id,"id"); pAttrib.add(id,"name"); 
    pAttrib.add(target,"target"); pAttrib.add(target,"ref"); 
    pAttrib.add(role,"role");
    pAttrib.put(cur);

    ParticleSet *qp = ptclPool->getParticleSet(target);
    if(qp == 0) {
      ERRORMSG("Wavefunction cannot be created because of missing particle set " << target)
      return false;
    }

    TrialWaveFunction *psi = getWaveFunction(id);
    if(psi) {
      WARNMSG("wavefunction with " << id << " is already created. Ignore the input")
      return true;
    }

    qp->setName(target);
    LOGMSG("Creating " << id << " wavefunction for " )
    LOGMSG(qp->getName() << " particleset")
    //Create a new TrialWaveFunction
    psi = new TrialWaveFunction;
    
    if(myPool.empty() || role == "primary") {
      primaryPsi=psi;
    }

    //Add to the pool
    myPool[id]=psi;

    cur = cur->children;
    while(cur != NULL) {
      string cname((const char*)(cur->name));
      if (cname == OrbitalBuilderBase::detset_tag) {
        string orbtype=(const char*)(xmlGetProp(cur, (const xmlChar *)"type"));
        LOGMSG("Slater-determinant terms using " << orbtype)
        if(orbtype == "MolecularOrbital") {
          MolecularOrbitalBuilder a(*qp,*psi,ptclPool->getPool());
          a.put(cur);
        } else {
          ERRORMSG(orbtype << " is disabled.")
          return false;
        }
      } else if (cname ==  OrbitalBuilderBase::jastrow_tag) {
        JastrowBuilder a(*qp,*psi,ptclPool->getPool());
        a.put(cur);
      }
      cur = cur->next;
    }
    return true;
  }

  bool WaveFunctionPool::put(std::istream& is) {
    return true;
  }

  bool WaveFunctionPool::get(std::ostream& os) const {
    return true;
  }

  void WaveFunctionPool::reset() {
 
  }
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
