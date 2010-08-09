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
#include "QMCWaveFunctions/ElectronGas/ElectronGasOrbitalBuilder.h"
#if defined(HAVE_EINSPLINE)
#include "QMCWaveFunctions/EinsplineSetBuilder.h"
#endif
#include "QMCWaveFunctions/OptimizableSPOBuilder.h"
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
//     delete_iter(basisBuilder.begin(),basisBuilder.end());
  }

  bool BasisSetFactory::put(xmlNodePtr cur) 
  {
    return true;
  }

  void BasisSetFactory::createBasisSet(xmlNodePtr cur,xmlNodePtr  rootNode) {

    ReportEngine PRE(ClassName,"createBasisSet");

    string sourceOpt("ion0");
    string typeOpt("");
    string name("");
    string keyOpt("NMO"); //gaussian Molecular Orbital
    string transformOpt("yes"); //numerical Molecular Orbital
    string cuspC("no");  // cusp correction
    string cuspInfo("");  // file with precalculated cusp correction info
    OhmmsAttributeSet aAttrib;
    aAttrib.add(sourceOpt,"source");
    aAttrib.add(cuspC,"cuspCorrection");
    aAttrib.add(typeOpt,"type");
    aAttrib.add(keyOpt,"keyword"); aAttrib.add(keyOpt,"key");
    aAttrib.add(name,"name");
    aAttrib.add(transformOpt,"transform");
    aAttrib.add(cuspInfo,"cuspInfo");
    if(rootNode != NULL)  aAttrib.put(rootNode); 
    
//     xmlNodePtr tc = cur->children; tc=tc->next;
//     while(tc != NULL) {
//     string cname;  getNodeName(cname,tc);
//     if (cname.find("asis")>1)
//     {
//       OhmmsAttributeSet bAttrib;
//       bAttrib.add(name,"name");
//       bAttrib.put(tc);
//       break;
//     }
//     tc=tc->next;
//     }

    BasisSetBuilder* bb=0;
    if(typeOpt.find("spline")<typeOpt.size())
    {
      name=typeOpt;
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
   
//      if(name=="") {
//        APP_ABORT("Missing basisset/@name.\n");
//      }

      if(transformOpt == "yes") { 
#if QMC_BUILD_LEVEL>2
        bb = new MolecularBasisBuilder<NGOBuilder>(targetPtcl,*ions,cuspC=="yes",cuspInfo);
#else
        bb = new MolecularBasisBuilder<NGOBuilder>(targetPtcl,*ions,false);
#endif
      } 
      else 
      {
#if QMC_BUILD_LEVEL>2
        if(cuspC == "yes")
           app_log() <<" ****** Cusp Correction algorithm is only implemented in combination with numerical radial orbitals. Use transform=yes to enable this option. \n";
#endif
        if(keyOpt == "GTO") 
          bb = new MolecularBasisBuilder<GTOBuilder>(targetPtcl,*ions);
        else if(keyOpt == "STO") 
          bb = new MolecularBasisBuilder<STOBuilder>(targetPtcl,*ions);
      }
    }
    else if (typeOpt == "linearopt") {
      //app_log()<<"Optimizable SPO set"<<endl;
      bb = new OptimizableSPOBuilder(targetPtcl,ptclPool,rootNode);
    }
    else if (typeOpt == "jellium") {
      app_log()<<"Electron gas SPO set"<<endl;
      bb = new ElectronGasBasisBuilder(targetPtcl,rootNode);
    }    

    PRE.flush();

    if(bb) 
    {
      bb->setReportLevel(ReportLevel);
      bb->initCommunicator(myComm);
      bb->put(cur);
      app_log()<<" Built basis "<< name<< endl;
//       basissets[name]=basisBuilder.size();
      basisBuilder[name]=(bb);
    } 
    else 
    {
      //fatal error
//       PRE.error("Failed to create a basis set.",true);
    }
  }

  SPOSetBase* BasisSetFactory::createSPOSet(xmlNodePtr cur)
  {
    string bname("");
    string bsname("");
    int bsize=basisBuilder.size(); 
    string sname(""); 
    OhmmsAttributeSet aAttrib; 
    aAttrib.add(bname,"basisset");
    aAttrib.add(bsname,"basis_sposet");
    aAttrib.add(sname,"name");
    aAttrib.put(cur); 
    
    string cname;
    xmlNodePtr tcur=cur->children;
    if (tcur!=NULL) getNodeName(cname,cur);
    
    if ( (basisBuilder.count(bname)==0 ) && (cname==basisset_tag)) 
      createBasisSet(tcur,cur);
    else if (basisBuilder.count(bsname))
    {
      createBasisSet(cur,cur);
      bname=sname;
    }
    else if (bname=="")
    {
      createBasisSet(cur,cur);
      bname=basisBuilder.rbegin()->first;
    }

    if(basisBuilder.size()) 
    {
      app_log()<<" Building SPOset "<<sname<<" with "<<bname<<" basis set."<<endl;
      return basisBuilder[bname]->createSPOSet(cur);
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
