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
#include "QMCWaveFunctions/JastrowBuilder.h"
#include "QMCWaveFunctions/MolecularOrbitals/MolecularOrbitalBuilder.h"
#include "QMCWaveFunctions/AtomicOrbitals/HeSTOClementiRottie.h"
#include "QMCWaveFunctions/ElectronGasOrbitalBuilder.h"
//#include "QMCWaveFunctions/NumericalOrbitalSetBuilder.h"
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
    if(psi == 0) {
      //Create a new TrialWaveFunction
      psi = new TrialWaveFunction;
      if(myPool.empty() || role == "primary") {
        primaryPsi=psi;
      }
      //Add to the pool
      myPool[id]=psi;
    } else {
      WARNMSG("wavefunction with " << id << " is already created. Add a new component.")
    }

    map<string,xmlNodePtr>::iterator wptr=m_wfsPtr.find(id);
    xmlNodePtr curWfsPtr=0;
    if(wptr == m_wfsPtr.end()) {
      curWfsPtr = xmlCopyNode(cur,2);
      m_wfsPtr[id]=curWfsPtr;
    } else {
      curWfsPtr=(*wptr).second;
    }

    cur = cur->children;
    while(cur != NULL) {
      string cname((const char*)(cur->name));
      if (cname == OrbitalBuilderBase::detset_tag) {
        string orbtype("MolecularOrbital");
        string nuclei("i");
        OhmmsAttributeSet oAttrib;
        oAttrib.add(orbtype,"type");
        oAttrib.add(nuclei,"source");
        oAttrib.put(cur);

        LOGMSG("Slater-determinant terms using " << orbtype)
        if(orbtype == "MolecularOrbital") {
          MolecularOrbitalBuilder a(*qp,*psi,ptclPool->getPool());
          a.put(cur);
        } else if(orbtype == "STO-He-Optimized") {
          ParticleSet* ion = ptclPool->getParticleSet(nuclei);
          if(ion) {
            HePresetHFBuilder a(*qp,*psi,*ion);
          }
        //} else if(orbtype == "NumericalOrbital") {
        //  NumericalOrbitalSetBuilder a(*qp,*psi,ptclPool->getPool());
        //  a.put(cur);
        } else if(orbtype == "electron-gas") {
          ElectronGasOrbitalBuilder a(*qp,*psi);
          a.put(cur);
        } else {
          ERRORMSG(orbtype << " is disabled.")
          return false;
        }
      } else if (cname ==  OrbitalBuilderBase::jastrow_tag) {
        xmlAddChild(curWfsPtr,xmlCopyNode(cur,1));
        JastrowBuilder a(*qp,*psi,ptclPool->getPool());
        a.put(cur);
      }
      cur = cur->next;
    }
    return true;
  }

  xmlNodePtr WaveFunctionPool::getWaveFunctionNode(const string& id) {

    if(m_wfsPtr.empty()) return 0;
    if(id == "null") {//return the first
      return (*(m_wfsPtr.begin())).second;
    } else {
      map<string,xmlNodePtr>::iterator wptr=m_wfsPtr.find(id);
      if(wptr == m_wfsPtr.end())
        return 0;
      else
        return (*wptr).second;
    }
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
