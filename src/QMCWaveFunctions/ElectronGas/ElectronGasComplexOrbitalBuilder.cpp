//////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Jeongnim Kim
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
#include "Utilities/OhmmsInfo.h"
#include "QMCWaveFunctions/ElectronGas/ElectronGasComplexOrbitalBuilder.h"
#include "QMCWaveFunctions/ElectronGas/HEGGrid.h"
#include "QMCWaveFunctions/Fermion/SlaterDet.h"
#include "OhmmsData/AttributeSet.h"

namespace qmcplusplus {

  /** constructor for EGOSet
   * @param norb number of orbitals for the EGOSet
   * @param k list of unique k points in Cartesian coordinate excluding gamma
   * @param k2 k2[i]=dot(k[i],k[i])
   */
  EGOSet::EGOSet(const vector<PosType>& k, const vector<RealType>& k2): K(k), mK2(k2) 
  {
    KptMax=k.size();
    Identity=true;
    OrbitalSetSize=k.size();
    BasisSetSize=k.size();
    className="EGOSet";
  }

  ElectronGasComplexOrbitalBuilder::ElectronGasComplexOrbitalBuilder(ParticleSet& els, 
      TrialWaveFunction& psi):
    OrbitalBuilderBase(els,psi)
  { 
  }   


  bool ElectronGasComplexOrbitalBuilder::put(xmlNodePtr cur){
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
    if(nc == 0) nc = egGrid.getShellIndex(nup);

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
    targetPsi.addOrbital(sdet,"SlaterDet");

    return true;
  }
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
