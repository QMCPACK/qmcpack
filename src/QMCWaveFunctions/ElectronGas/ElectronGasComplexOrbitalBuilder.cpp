//////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Jeongnim Kim
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
#include "Utilities/OhmmsInfo.h"
#include "QMCWaveFunctions/ElectronGas/ElectronGasComplexOrbitalBuilder.h"
#include "QMCWaveFunctions/ElectronGas/HEGGrid.h"
#include "QMCWaveFunctions/SlaterDeterminant.h"

namespace qmcplusplus {

  /** constructor for EGOSet
   * @param norb number of orbitals for the EGOSet
   * @param k list of unique k points in Cartesian coordinate excluding gamma
   * @param k2 k2[i]=dot(k[i],k[i])
   */
  ElectronGasComplexOrbitalBuilder::EGOSet::EGOSet(const vector<PosType>& k, 
      const vector<RealType>& k2): K(k), mK2(k2) {
    KptMax=k.size();
  }

  ElectronGasComplexOrbitalBuilder::ElectronGasComplexOrbitalBuilder(ParticleSet& els, 
      TrialWaveFunction& psi):
    OrbitalBuilderBase(els,psi)
  { 
  }   


  bool ElectronGasComplexOrbitalBuilder::put(xmlNodePtr cur){

    //can use generic AttributeSet but leave it for now
    int nc=0;
    const xmlChar* nc_ptr=xmlGetProp(cur,(const xmlChar*)"shell");
    if(nc_ptr) {
      nc = atoi((const char*)nc_ptr);
    }

    PosType twistAngle;
    nc_ptr=xmlGetProp(cur,(const xmlChar*)"twist");
    if(nc_ptr) {
      //putAttribute(shift,nc_ptr);
      std::istringstream stream((const char*)nc_ptr);
      stream >> twistAngle;
    }

    typedef DiracDeterminant<EGOSet>  Det_t;
    typedef SlaterDeterminant<EGOSet> SlaterDeterminant_t;

    int nat=targetPtcl.getTotalNum();
    int nup=nat/2;

    HEGGrid<RealType,OHMMS_DIM> egGrid(targetPtcl.Lattice);
    if(nc == 0) {
      nc = egGrid.getNC(nup);
    }

    //number of kpoints in a half sphere at zero twist
    int nkpts=(nup-1)/2;

    //create k-points ordered wrt the k^2
    egGrid.createGrid(nc,nkpts);
    egGrid.createGrid(twistAngle);

    //create a E(lectron)G(as)O(rbital)Set
    EGOSet* psi=new EGOSet(egGrid.kpt,egGrid.mk2); 

    //create up determinant
    Det_t *updet = new Det_t(*psi,0);
    updet->set(0,nup);

    //create down determinant
    Det_t *downdet = new Det_t(*psi,nup);
    downdet->set(nup,nup);

    //create a Slater determinant
    SlaterDeterminant_t *sdet  = new SlaterDeterminant_t;
    sdet->add(updet);
    sdet->add(downdet);

    //add a DummyBasisSet
    sdet->setBasisSet(new DummyBasisSet);

    //add Slater determinant to targetPsi
    targetPsi.addOrbital(sdet);

    return true;
  }
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
