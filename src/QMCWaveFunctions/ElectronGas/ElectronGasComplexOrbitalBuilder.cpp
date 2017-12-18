//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#include "QMCWaveFunctions/ElectronGas/ElectronGasComplexOrbitalBuilder.h"
#include "QMCWaveFunctions/Fermion/SlaterDet.h"
#if QMC_BUILD_LEVEL>2
#include "QMCWaveFunctions/Fermion/BackflowTransformation.h"
#endif
#include "OhmmsData/AttributeSet.h"

namespace qmcplusplus
{

/** constructor for EGOSet
 * @param norb number of orbitals for the EGOSet
 * @param k list of unique k points in Cartesian coordinate excluding gamma
 * @param k2 k2[i]=dot(k[i],k[i])
 */
EGOSet::EGOSet(const std::vector<PosType>& k, const std::vector<RealType>& k2): K(k), mK2(k2)
{
  KptMax=k.size();
  Identity=true;
  OrbitalSetSize=k.size();
  BasisSetSize=k.size();
  className="EGOSet";
  //assign_energies();
}

EGOSet::EGOSet(const std::vector<PosType>& k, const std::vector<RealType>& k2, const std::vector<int>& d)
  : K(k), mK2(k2)
{
  KptMax=k.size();
  Identity=true;
  OrbitalSetSize=k.size();
  BasisSetSize=k.size();
  className="EGOSet";
  //assign_energies();
  //assign_degeneracies(d);
}

ElectronGasComplexOrbitalBuilder::ElectronGasComplexOrbitalBuilder(ParticleSet& els,
    TrialWaveFunction& psi):
  OrbitalBuilderBase(els,psi)
{
}


bool ElectronGasComplexOrbitalBuilder::put(xmlNodePtr cur)
{
  int nc=0;
  PosType twist(0.0);
  OhmmsAttributeSet aAttrib;
  aAttrib.add(nc,"shell");
  aAttrib.add(twist,"twist");
  aAttrib.put(cur);
  //typedef DiracDeterminant<EGOSet>  Det_t;
  //typedef SlaterDeterminant<EGOSet> SlaterDeterminant_t;
  typedef DiracDeterminantBase  Det_t;
  typedef SlaterDet SlaterDeterminant_t;
  int nat=targetPtcl.getTotalNum();
  int nup=nat/2;
  HEGGrid<RealType,OHMMS_DIM> egGrid(targetPtcl.Lattice);
  if(nc == 0)
    nc = egGrid.getShellIndex(nup);
  egGrid.createGrid(nc,nup,twist);
  targetPtcl.setTwist(twist);
  //create a E(lectron)G(as)O(rbital)Set
  EGOSet* psiu=new EGOSet(egGrid.kpt,egGrid.mk2);
  EGOSet* psid=new EGOSet(egGrid.kpt,egGrid.mk2);
  //create up determinant
  Det_t *updet = new Det_t(psiu);
  updet->set(0,nup);
  //create down determinant
  Det_t *downdet = new Det_t(psid);
  downdet->set(nup,nup);
  //create a Slater determinant
  //SlaterDeterminant_t *sdet  = new SlaterDeterminant_t;
  SlaterDet *sdet  = new SlaterDet(targetPtcl);
  sdet->add(psiu,"u");
  sdet->add(psid,"d");
  sdet->add(updet,0);
  sdet->add(downdet,1);
  //add Slater determinant to targetPsi
  targetPsi.addOrbital(sdet,"SlaterDet",true);
  return true;
}

ElectronGasBasisBuilder::ElectronGasBasisBuilder(ParticleSet& p, xmlNodePtr cur)
  :egGrid(p.Lattice),unique_twist(-1.0),has_twist(false)
{
}

bool ElectronGasBasisBuilder::put(xmlNodePtr cur)
{
  OhmmsAttributeSet aAttrib;
  aAttrib.add(unique_twist,"twist");
  aAttrib.put(cur);

  has_twist = true;
  for(int d=0;d<OHMMS_DIM;++d)
    has_twist &= (unique_twist[d]+1.0)>1e-6;

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
  if(has_twist)
    twist = unique_twist;
  if(ns>0)
    nc = egGrid.getShellIndex(ns);
  if (nc<0)
  {
    app_error() << "  HEG Invalid Shell." << std::endl;
    APP_ABORT("ElectronGasBasisBuilder::put");
  }
  egGrid.createGrid(nc,ns,twist);
  EGOSet* spo = new EGOSet(egGrid.kpt,egGrid.mk2,egGrid.deg);
  return spo;
}


SPOSetBase* ElectronGasBasisBuilder::createSPOSetFromIndices(indices_t& indices)
{
  egGrid.createGrid(indices);
  EGOSet* spo = new EGOSet(egGrid.kpt,egGrid.mk2,egGrid.deg);
  return spo;
}






}
