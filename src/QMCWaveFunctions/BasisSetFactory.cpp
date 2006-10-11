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
#include "QMCWaveFunctions/BasisSetFactory.h"
//#include "QMCWaveFunctions/MolecularOrbitals/Any2GridBuilder.h"
#include "QMCWaveFunctions/MolecularOrbitals/NGOBuilder.h"
#include "QMCWaveFunctions/MolecularOrbitals/GTOBuilder.h"
#include "QMCWaveFunctions/MolecularOrbitals/STOBuilder.h"
#include "QMCWaveFunctions/MolecularOrbitals/MolecularBasisBuilder.h"
#include "QMCWaveFunctions/SplineSetBuilder.h"
#include "Message/Communicate.h"
#include "OhmmsData/AttributeSet.h"

namespace qmcplusplus {

  /** constructor
   * \param els reference to the electrons
   * \param psi reference to the wavefunction
   * \param ions reference to the ions
   */
  BasisSetFactory::BasisSetFactory(ParticleSet& els, TrialWaveFunction& psi, PtclPoolType& psets):
  OrbitalBuilderBase(els,psi), ptclPool(psets){

  }

  bool BasisSetFactory::put(xmlNodePtr cur) {
    return true;
  }

  void BasisSetFactory::createBasisSet(xmlNodePtr cur, xmlNodePtr rootNode) {

    string sourceOpt("ion0");
    string typeOpt("MolecularOrbital");
    string keyOpt("NMO"); //numerical Molecular Orbital
    string transformOpt("yes"); //numerical Molecular Orbital
    OhmmsAttributeSet aAttrib;
    aAttrib.add(sourceOpt,"source");
    aAttrib.add(typeOpt,"type");
    aAttrib.add(keyOpt,"keyword"); aAttrib.add(keyOpt,"key");
    aAttrib.add(transformOpt,"transform");
    if(rootNode != NULL)  {
      aAttrib.put(rootNode);
    }

    BasisSetBuilder* bb=0;
    if(typeOpt == "spline") {
      app_log() << "  SplineSetBuilder: spline on 3D TriCubicGrid " << endl;
      bb = new SplineSetBuilder(targetPtcl,ptclPool);
    } else { 
    //if(typeOpt == "MolecularOrbital") {
      ParticleSet* ions=0;

      //initialize with the source tag
      PtclPoolType::iterator pit(ptclPool.find(sourceOpt));
      if(pit == ptclPool.end()) {
        app_error() << "Molecular orbital cannot be created.\n" 
          << "Missing/incorrect source attribute. Abort at BasisSetFactory::createBasisSet." << endl;
        OHMMS::Controller->abort();
      } else {
        app_log() << "  Molecular orbital with " << sourceOpt;
        ions=(*pit).second; 
      }

      if(transformOpt == "yes") {
        app_log() << " by numerical radial functors." << endl;
        bb = new MolecularBasisBuilder<NGOBuilder>(targetPtcl,*ions);
        //bb = new MolecularBasisBuilder<Any2GridBuilder>(targetPtcl,*ions);
      } else {
        if(keyOpt == "GTO") {
          app_log() << " by analytic GTO functors." << endl;
          bb = new MolecularBasisBuilder<GTOBuilder>(targetPtcl,*ions);
        } else if(keyOpt == "STO") {
          app_log() << " by analytic STO functors." << endl;
          bb = new MolecularBasisBuilder<STOBuilder>(targetPtcl,*ions);
        }
      }
    }

    if(bb) {
      bb->put(cur);
      basisBuilder.push_back(bb);
    } else {
      app_log() << endl;
      app_error() << "  Failed to create a basis set. Stop at BasisSetFactory::createBasisSet" << endl;
      OHMMS::Controller->abort();
    }
  }

  SPOSetBase* BasisSetFactory::createSPOSet(xmlNodePtr cur) {
    if(basisBuilder.size()) {
      return basisBuilder.back()->createSPOSet(cur);
    } else {
      app_error() << "  Failed to create a SPOSet. Stop at BasisSetFactory::createBasisSet" << endl;
      OHMMS::Controller->abort();
      return 0;
    }
  }
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
