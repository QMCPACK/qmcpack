//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#include "QMCWaveFunctions/ElectronGas/ElectronGasOrbitalBuilder.h"
#include "QMCWaveFunctions/Fermion/SlaterDet.h"
#include "QMCWaveFunctions/Fermion/RNDiracDeterminantBase.h"
#include "QMCWaveFunctions/Fermion/RNDiracDeterminantBaseAlternate.h"
#include "OhmmsData/AttributeSet.h"
#if QMC_BUILD_LEVEL>2
#include "QMCWaveFunctions/Fermion/BackflowBuilder.h"
#include "QMCWaveFunctions/Fermion/BackflowTransformation.h"
#include "QMCWaveFunctions/Fermion/SlaterDetWithBackflow.h"
#include "QMCWaveFunctions/Fermion/DiracDeterminantWithBackflow.h"
#endif

namespace qmcplusplus
{

/** constructor for EGOSet
 * @param norb number of orbitals for the EGOSet
 * @param k list of unique k points in Cartesian coordinate excluding gamma
 * @param k2 k2[i]=dot(k[i],k[i])
 */
RealEGOSet::RealEGOSet(const std::vector<PosType>& k, const std::vector<RealType>& k2): K(k),mK2(k2)
{
  KptMax=k.size();
  Identity=true;
  OrbitalSetSize=2*k.size()+1;
  BasisSetSize=2*k.size()+1;
  className="EGOSet";
}

ElectronGasOrbitalBuilder::ElectronGasOrbitalBuilder(ParticleSet& els, TrialWaveFunction& psi):
#if QMC_BUILD_LEVEL>2
  OrbitalBuilderBase(els,psi),UseBackflow(false),BFTrans(0)
#else
  OrbitalBuilderBase(els,psi),UseBackflow(false)
#endif
{
}

bool ElectronGasOrbitalBuilder::put(xmlNodePtr cur)
{
  int nc(0), nc2(-2);
  ValueType pol(0);
  ValueType bosonic_eps(-999999);
  ValueType rntype(0);
  PosType twist(0.0);
  OhmmsAttributeSet aAttrib;
  aAttrib.add(nc,"shell");
  aAttrib.add(nc2,"shell2");
  aAttrib.add(bosonic_eps,"eps");
  aAttrib.add(rntype,"primary");
  aAttrib.add(twist,"twist");
  aAttrib.put(cur);
  if (nc2==-2)
    nc2=nc;
#if QMC_BUILD_LEVEL>2
  xmlNodePtr curRoot=cur;
  xmlNodePtr BFNode(NULL);
  std::string cname;
  cur = curRoot->children;
  while (cur != NULL)//check the basis set
  {
    getNodeName(cname,cur);
    if(cname == backflow_tag)
    {
// FIX FIX FIX !!!
      UseBackflow=true;
      if(BFNode == NULL)
        BFNode=cur;
    }
    cur = cur->next;
  }
#endif
  typedef SlaterDet::Determinant_t Det_t;
  typedef SlaterDet SlaterDeterminant_t;
  HEGGrid<RealType,OHMMS_DIM> egGrid(targetPtcl.Lattice);
  HEGGrid<RealType,OHMMS_DIM> egGrid2(targetPtcl.Lattice);
  int nat=targetPtcl.getTotalNum();
  if (nc == 0)
    nc = nc2 = egGrid.getShellIndex(nat/2);
  int nup= egGrid.getNumberOfKpoints(nc);
  int ndn(0);
  if (nc2>-1)
    ndn = egGrid.getNumberOfKpoints(nc2);
  if (nc<0)
  {
    app_error() << "  HEG Invalid Shell." << std::endl;
    APP_ABORT("ElectronGasOrbitalBuilder::put");
  }
  if (nat!=(nup+ndn))
  {
    app_error() << "  The number of particles "<<nup<<"/"<<ndn<< " does not match to the shell." << std::endl;
    app_error() << "  Suggested values for the number of particles " << std::endl;
    app_error() << "   " << 2*egGrid.getNumberOfKpoints(nc) << " for shell "<< nc << std::endl;
    app_error() << "   " << 2*egGrid.getNumberOfKpoints(nc-1) << " for shell "<< nc-1 << std::endl;
    APP_ABORT("ElectronGasOrbitalBuilder::put");
    return false;
  }
  int nkpts=(nup-1)/2;
  int nkpts2=(ndn-1)/2;
  RealEGOSet* psiu;
  RealEGOSet* psid;
  if(nup==ndn)
  {
    //create a E(lectron)G(as)O(rbital)Set
    egGrid.createGrid(nc,nkpts);
    psiu=new RealEGOSet(egGrid.kpt,egGrid.mk2);
    psid=new RealEGOSet(egGrid.kpt,egGrid.mk2);
  }
  else
    if(ndn>0)
    {
      //create a E(lectron)G(as)O(rbital)Set
      egGrid.createGrid(nc,nkpts);
      egGrid2.createGrid(nc2,nkpts2);
      psiu=new RealEGOSet(egGrid.kpt,egGrid.mk2);
      psid=new RealEGOSet(egGrid.kpt,egGrid.mk2);
    }
    else
    {
      //create a E(lectron)G(as)O(rbital)Set
      egGrid.createGrid(nc,nkpts);
      psiu=new RealEGOSet(egGrid.kpt,egGrid.mk2);
    }
  //create a Slater determinant
  SlaterDeterminant_t *sdet;
#if QMC_BUILD_LEVEL>2
  if(UseBackflow)
    sdet = new SlaterDetWithBackflow(targetPtcl,BFTrans);
  else
#endif
    sdet  = new SlaterDeterminant_t(targetPtcl);
  //add SPOSets
  sdet->add(psiu,"u");
  if(ndn>0)
    sdet->add(psid,"d");
  if (rntype>0)
  {
#if QMC_BUILD_LEVEL>2
    if(UseBackflow)
      APP_ABORT("RN with Backflow not implemented. \n");
#endif
    if (rntype==1)
      app_log()<<" Using determinant with eps="<<bosonic_eps<< std::endl;
    else
      app_log()<<" Using alternate determinant with eps="<<bosonic_eps<< std::endl;
    //create up determinant
    Det_t *updet;
    if (rntype==1)
      updet = new RNDiracDeterminantBase(psiu);
    else
      updet = new RNDiracDeterminantBaseAlternate(psiu);
    updet->setLogEpsilon(bosonic_eps);
    updet->set(0,nup);
    sdet->add(updet,0);
    if(ndn>0)
    {
      //create down determinant
      Det_t *downdet;
      if (rntype==1)
        downdet = new RNDiracDeterminantBase(psiu);
      else
        downdet = new RNDiracDeterminantBaseAlternate(psiu);
      downdet->set(nup,ndn);
      downdet->setLogEpsilon(bosonic_eps);
      sdet->add(downdet,1);
    }
  }
  else
  {
    Det_t *updet, *downdet;
#if QMC_BUILD_LEVEL>2
    if(UseBackflow)
    {
      app_log() <<"Creating Backflow transformation in ElectronGasOrbitalBuilder::put(xmlNodePtr cur).\n";
      //create up determinant
      updet = new DiracDeterminantWithBackflow(targetPtcl,psiu,BFTrans,0);
      updet->set(0,nup);
      if(ndn>0)
      {
        //create down determinant
        downdet = new DiracDeterminantWithBackflow(targetPtcl,psid,BFTrans,nup);
        downdet->set(nup,ndn);
      }
      PtclPoolType dummy;
      BackflowBuilder* bfbuilder = new BackflowBuilder(targetPtcl,dummy,targetPsi);
      bfbuilder->put(BFNode);
      BFTrans = bfbuilder->getBFTrans();
      sdet->add(updet,0);
      if(ndn>0)
        sdet->add(downdet,1);
      sdet->setBF(BFTrans);
      if(BFTrans->isOptimizable())
        sdet->Optimizable = true;
      sdet->resetTargetParticleSet(targetPtcl);
    }
    else
#endif
    {
      //create up determinant
      updet = new Det_t(psiu);
      updet->set(0,nup);
      if(ndn>0)
      {
        //create down determinant
        downdet = new Det_t(psid);
        downdet->set(nup,ndn);
      }
      sdet->add(updet,0);
      if(ndn>0)
        sdet->add(downdet,1);
    }
  }
//#if QMC_BUILD_LEVEL>2
//    // change DistanceTables if using backflow
//    if(UseBackflow)   {
//       sdet->resetTargetParticleSet(targetPtcl);
//    }
//#endif
  //add Slater determinant to targetPsi
  targetPsi.addOrbital(sdet,"SlaterDet",true);
  return true;
}

ElectronGasBasisBuilder::ElectronGasBasisBuilder(ParticleSet& p, xmlNodePtr cur)
  :egGrid(p.Lattice)
{
}

bool ElectronGasBasisBuilder::put(xmlNodePtr cur)
{
  return true;
}
  
SPOSetBase* ElectronGasBasisBuilder::createSPOSetFromXML(xmlNodePtr cur)
{
  app_log() << "ElectronGasBasisBuilder::createSPOSet " << std::endl;
  int nc=0;
  int ns=0;
  PosType twist(0.0);
  std::string spo_name("heg");
  OhmmsAttributeSet aAttrib;
  aAttrib.add(ns,"size");
  aAttrib.add(twist,"twist");
  aAttrib.add(spo_name,"name");
  aAttrib.add(spo_name,"id");
  aAttrib.put(cur);
  if (ns>0)
    nc=egGrid.getShellFromStates(ns);
  if (nc<0)
  {
    app_error() << "  HEG Invalid Shell." << std::endl;
    APP_ABORT("ElectronGasOrbitalBuilder::put");
  }
  egGrid.createGrid(nc,(ns-1)/2);
  return new RealEGOSet(egGrid.kpt,egGrid.mk2);
}

}
