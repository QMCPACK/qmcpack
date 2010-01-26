//////////////////////////////////////////////////////////////////
// (c) Copyright 2006-  by Jeongnim Kim
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
#include "QMCWaveFunctions/BasisSetFactory.h"
#if OHMMS_DIM == 3
#include "QMCWaveFunctions/MolecularOrbitals/NGOBuilder.h"
#include "QMCWaveFunctions/MolecularOrbitals/GTOBuilder.h"
#include "QMCWaveFunctions/MolecularOrbitals/STOBuilder.h"
#include "QMCWaveFunctions/MolecularOrbitals/MolecularBasisBuilder.h"
#if defined(HAVE_EINSPLINE)
#include "QMCWaveFunctions/EinsplineSetBuilder.h"
#endif
#if QMC_BUILD_LEVEL>1
#include "QMCWaveFunctions/TricubicBsplineSetBuilder.h"
#endif
#endif
#include "Utilities/ProgressReportEngine.h"
#include "Utilities/IteratorUtility.h"
#include "OhmmsData/AttributeSet.h"

namespace qmcplusplus {

  /** constructor
   * \param els reference to the electrons
   * \param psi reference to the wavefunction
   * \param ions reference to the ions
   */
  BasisSetFactory::BasisSetFactory(ParticleSet& els, TrialWaveFunction& psi, PtclPoolType& psets):
  OrbitalBuilderBase(els,psi), ptclPool(psets)
  {
    ClassName="BasisSetFactory";
  }

  BasisSetFactory::~BasisSetFactory()
  {
    DEBUG_MEMORY("BasisSetFactory::~BasisSetFactory");
    delete_iter(basisBuilder.begin(),basisBuilder.end());
  }

  bool BasisSetFactory::put(xmlNodePtr cur) 
  {
    return true;
  }

  void BasisSetFactory::createBasisSet(xmlNodePtr cur, xmlNodePtr rootNode) {

    ReportEngine PRE(ClassName,"createBasisSet");

    string sourceOpt("ion0");
    string typeOpt("MolecularOrbital");
    string keyOpt("NMO"); //numerical Molecular Orbital
    string transformOpt("yes"); //numerical Molecular Orbital
    OhmmsAttributeSet aAttrib;
    aAttrib.add(sourceOpt,"source");
    aAttrib.add(typeOpt,"type");
    aAttrib.add(keyOpt,"keyword"); aAttrib.add(keyOpt,"key");
    aAttrib.add(transformOpt,"transform");
    if(rootNode != NULL)  aAttrib.put(rootNode);

    BasisSetBuilder* bb=0;
    if(typeOpt.find("spline")<typeOpt.size())
    {
#if defined(HAVE_EINSPLINE)
      PRE << "EinsplineSetBuilder:  using libeinspline for B-spline orbitals.\n";
      bb = new EinsplineSetBuilder(targetPtcl,ptclPool,rootNode);
#else
      PRE.error("Einspline is missing for B-spline orbitals",true);
      //PRE << "TricubicBsplineSetBuilder: b-spline on 3D TriCubicGrid.\n";
      //bb = new TricubicBsplineSetBuilder(targetPtcl,ptclPool,rootNode);
#endif
    }
    else if(typeOpt == "MolecularOrbital" || typeOpt == "MO") 
    {
      ParticleSet* ions=0;

      //do not use box to check the boundary conditions
      if(targetPtcl.Lattice.SuperCellEnum==SUPERCELL_OPEN) targetPtcl.setBoundBox(false);

      //initialize with the source tag
      PtclPoolType::iterator pit(ptclPool.find(sourceOpt));
      if(pit == ptclPool.end()) 
        PRE.error("Missing basisset/@source.",true);
      else 
        ions=(*pit).second; 

      if(transformOpt == "yes") 
        bb = new MolecularBasisBuilder<NGOBuilder>(targetPtcl,*ions);
      else 
      {
        if(keyOpt == "GTO") 
          bb = new MolecularBasisBuilder<GTOBuilder>(targetPtcl,*ions);
        else if(keyOpt == "STO") 
          bb = new MolecularBasisBuilder<STOBuilder>(targetPtcl,*ions);
      }
    }

    PRE.flush();

    if(bb) 
    {
      bb->setReportLevel(ReportLevel);
      bb->initCommunicator(myComm);
      bb->put(cur);
      basisBuilder.push_back(bb);
    } 
    else 
    {
      //fatal error
      PRE.error("Failed to create a basis set.",true);
    }
  }

  SPOSetBase* BasisSetFactory::createSPOSet(xmlNodePtr cur) 
  {
    if(basisBuilder.size()) 
    {
      return basisBuilder.back()->createSPOSet(cur);
    } 
    else 
    {
      APP_ABORT("BasisSetFactory::createSPOSet Failed to create a SPOSet. basisBuilder is empty.");
      return 0;
    }
  }
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
