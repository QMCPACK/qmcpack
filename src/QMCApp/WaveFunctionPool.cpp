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

    string id("primary"), target("e");
    OhmmsAttributeSet pAttrib;
    pAttrib.add(id,"id"); pAttrib.add(id,"name"); 
    pAttrib.add(target,"target"); pAttrib.add(target,"ref"); 

    ParticleSet *p = ptclPool->getParticleSet(target);
    if(p == 0) {
      ERRORMSG("Wavefunction cannot be created because of missing particle set " << target)
      return false;
    }

    TrialWaveFunction *psi = getWaveFunction(id);

    if(psi) {
      WARNMSG("wavefunction with " << id << " is already created. Ignore the input")
      return true;
    }

    psi = new TrialWaveFunction;

    //add particles associated with this wave function
    cur = cur->children;
    while(cur != NULL) {
      string cname((const char*)(cur->name));
      if (cname == OrbitalBuilderBase::detset_tag) {
        string orbtype=(const char*)(xmlGetProp(cur, (const xmlChar *)"type"));
        LOGMSG("Slater-determinant terms using " << orbtype)
        if(orbtype == "MolecularOrbital") {
          MolecularOrbitalBuilder a(*p,*psi,ptclPool->getPool());
          a.put(cur);
        }
      } else if (cname ==  OrbitalBuilderBase::jastrow_tag) {
        JastrowBuilder a(*p,*psi,ptclPool->getPool());
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
