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
#include "QMCWaveFunctions/ElectronGas/ElectronGasOrbitalBuilder.h"
#include "QMCWaveFunctions/HarmonicOscillator/SHOSetBuilder.h"
#if OHMMS_DIM == 3
#if !defined(QMC_COMPLEX)
#include "QMCWaveFunctions/MolecularOrbitals/NGOBuilder.h"
#include "QMCWaveFunctions/MolecularOrbitals/GTOBuilder.h"
#include "QMCWaveFunctions/MolecularOrbitals/STOBuilder.h"
#include "QMCWaveFunctions/MolecularOrbitals/MolecularBasisBuilder.h"
#endif

#if defined(HAVE_EINSPLINE)
#include "QMCWaveFunctions/EinsplineSetBuilder.h"
#endif
#endif
#include "QMCWaveFunctions/OptimizableSPOBuilder.h"
#include "QMCWaveFunctions/AFMSPOBuilder.h"
#include "Utilities/ProgressReportEngine.h"
#include "Utilities/IteratorUtility.h"
#include "OhmmsData/AttributeSet.h"

namespace qmcplusplus
{


  SPOSetBase* get_sposet(const string& name)
  {
    int nfound = 0;
    SPOSetBase* spo = 0;
    map<string,BasisSetBuilder*>::iterator it;
    for(it=basis_builders.begin();it!=basis_builders.end();++it)
    {
      vector<SPOSetBase*>& sposets = it->second->sposets;
      for(int i=0;i<sposets.size();++i)
      {
        SPOSetBase* sposet = sposets[i];
        if(sposet->objectName==name)
        {
          spo = sposet;
          nfound++;
        }
      }
    }
    if(nfound>1)
    {
      write_basis_builders();
      APP_ABORT("get_sposet: requested sposet "+name+" is not unique");
    }
    else if(spo==NULL)
    {
      write_basis_builders();
      APP_ABORT("get_sposet: requested sposet "+name+" does not exist");
    }
    return spo;
  }


  void write_basis_builders(const string& pad)
  {
    string pad2 = pad+"  ";
    map<string,BasisSetBuilder*>::iterator it;
    for(it=basis_builders.begin();it!=basis_builders.end();++it)
    {
      const string& type = it->first;
      vector<SPOSetBase*>& sposets = it->second->sposets;
      app_log()<<pad<<"sposets for BasisSetBuilder of type "<<type<<endl;
      for(int i=0;i<sposets.size();++i)
      {
        app_log()<<pad2<<"sposet "<<sposets[i]->objectName<<endl;
      }
    }
  }



/** constructor
 * \param els reference to the electrons
 * \param psi reference to the wavefunction
 * \param ions reference to the ions
 */
BasisSetFactory::BasisSetFactory(ParticleSet& els, TrialWaveFunction& psi, PtclPoolType& psets):
  OrbitalBuilderBase(els,psi), ptclPool(psets), last_builder(0)
{
  ClassName="BasisSetFactory";
}

BasisSetFactory::~BasisSetFactory()
{
  DEBUG_MEMORY("BasisSetFactory::~BasisSetFactory");
}

bool BasisSetFactory::put(xmlNodePtr cur)
{
  return true;
}

BasisSetBuilder* BasisSetFactory::createBasisSet(xmlNodePtr cur,xmlNodePtr  rootNode)
{
  ReportEngine PRE(ClassName,"createBasisSet");
  string sourceOpt("ion0");
  string type("");
  string name("");
  string keyOpt("NMO"); //gaussian Molecular Orbital
  string transformOpt("yes"); //numerical Molecular Orbital
  string cuspC("no");  // cusp correction
  string cuspInfo("");  // file with precalculated cusp correction info
  OhmmsAttributeSet aAttrib;
  aAttrib.add(sourceOpt,"source");
  aAttrib.add(cuspC,"cuspCorrection");
  aAttrib.add(type,"type");
  aAttrib.add(keyOpt,"keyword");
  aAttrib.add(keyOpt,"key");
  aAttrib.add(name,"name");
  aAttrib.add(transformOpt,"transform");
  aAttrib.add(cuspInfo,"cuspInfo");
  if(rootNode != NULL)
    aAttrib.put(rootNode);

  app_log()<<"  basis builder: type = "<<type<<endl;

  tolower(type);

  app_log()<<"  basis builder: type = "<<type<<endl;

  BasisSetBuilder* bb=0;
  if (type == "jellium" || type == "heg")
  {
    app_log()<<"Electron gas SPO set"<<endl;
    bb = new ElectronGasBasisBuilder(targetPtcl,rootNode);
  }
  else if (type == "sho")
  {
    app_log()<<"Harmonic Oscillator SPO set"<<endl;
    bb = new SHOSetBuilder(targetPtcl);
  }
  else if (type == "linearopt")
  {
    //app_log()<<"Optimizable SPO set"<<endl;
    bb = new OptimizableSPOBuilder(targetPtcl,ptclPool,rootNode);
  }
  else if (type == "afm")
  {
    //       app_log()<<"AFM SPO set"<<endl;
    bb = new AFMSPOBuilder(targetPtcl,ptclPool,rootNode);
  }
#if OHMMS_DIM ==3
  else if(type.find("spline")<type.size())
  {
    name=type;
#if defined(HAVE_EINSPLINE)
    PRE << "EinsplineSetBuilder:  using libeinspline for B-spline orbitals.\n";
    bb = new EinsplineSetBuilder(targetPtcl,ptclPool,rootNode);
#else
    PRE.error("Einspline is missing for B-spline orbitals",true);
#endif
  }
#if !defined(QMC_COMPLEX)
  else if(type == "molecularorbital" || type == "mo")
  {
    ParticleSet* ions=0;
    //do not use box to check the boundary conditions
    if(targetPtcl.Lattice.SuperCellEnum==SUPERCELL_OPEN)
      targetPtcl.setBoundBox(false);
    //initialize with the source tag
    PtclPoolType::iterator pit(ptclPool.find(sourceOpt));
    if(pit == ptclPool.end())
      PRE.error("Missing basisset/@source.",true);
    else
      ions=(*pit).second;
    if(transformOpt == "yes")
    {
      app_log() << "Using MolecularBasisBuilder<NGOBuilder>" << endl;
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
  else 
  {
    APP_ABORT("BasisSetFactory::createBasisSet\n basis builder type="+type+" is unknown");
  }
#endif //!QMC_COMPLEX
#endif  //OHMMS_DIM==3
  PRE.flush();
  if(bb)
  {
    bb->setReportLevel(ReportLevel);
    bb->initCommunicator(myComm);
    bb->put(cur);
    app_log()<<" Built BasisSetBuilder of type "<< type << endl;
    basis_builders[type] = bb;
  }
  else
    APP_ABORT("BasisSetFactory::createBasisSet\n  BasisSetBuilder creation failed.");

  last_builder = bb;

  return bb;
}


SPOSetBase* BasisSetFactory::createSPOSet(xmlNodePtr cur)
{
  string bname("");
  string bsname("");
  string sname("");
  string type("");
  OhmmsAttributeSet aAttrib;
  aAttrib.add(bname,"basisset");
  aAttrib.add(bsname,"basis_sposet");
  aAttrib.add(sname,"name");
  aAttrib.add(type,"type");
  //aAttrib.put(rcur);
  aAttrib.put(cur);
  tolower(type);


  BasisSetBuilder* bb;

  if(bname=="")
    bname=type;
  if(type=="")
    bb = last_builder;
  else if(basis_builders.find(type)!=basis_builders.end())
    bb = basis_builders[type];
  else
  {
    string cname("");
    xmlNodePtr tcur=cur->children;
    if(tcur!=NULL)
      getNodeName(cname,cur);
    if(cname==basisset_tag)
      bb = createBasisSet(tcur,cur);
    else
      bb = createBasisSet(cur,cur);
  }
  if(bb)
  {
    app_log()<<" Building SPOset "<<sname<<" with "<<bname<<" basis set."<<endl;
    return bb->createSPOSet(cur);
  }
  else
  {
    APP_ABORT("BasisSetFactory::createSPOSet Failed to create a SPOSet. basisBuilder is empty.");
    return 0;
  }
}


void BasisSetFactory::build_sposet_collection(xmlNodePtr cur)
{
  xmlNodePtr parent = cur;
  string type("");
  OhmmsAttributeSet attrib;
  attrib.add(type,"type");
  attrib.put(cur);
  tolower(type);

  app_log()<<"building sposet collection of type "<<type<<endl;

  BasisSetBuilder* bb = createBasisSet(cur,cur);
  xmlNodePtr element = parent->children;
  int nsposets = 0;
  while(element!=NULL)
  {
    string cname((const char*)(element->name));
    if(cname=="sposet")
    {
      string name("");
      OhmmsAttributeSet attrib;
      attrib.add(name,"name");
      attrib.put(element);

      app_log()<<"  Building SPOSet "<<name<<" with "<<type<<" BasisSetBuilder"<<endl;
      SPOSetBase* spo = bb->createSPOSet(element);
      spo->objectName = name;
      nsposets++;
    }
    element = element->next;
  }
  if(nsposets==0)
    APP_ABORT("BasisSetFactory::build_sposet_collection  no <sposet/> elements found");
}


}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/
