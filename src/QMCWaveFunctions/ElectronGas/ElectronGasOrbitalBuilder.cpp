/////////////////////////////////////////////////////////////////
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
#include "QMCWaveFunctions/ElectronGas/ElectronGasOrbitalBuilder.h"
#include "QMCWaveFunctions/SlaterDeterminant.h"
#include "QMCWaveFunctions/ElectronGas/HEGGrid.h"
#include "OhmmsData/AttributeSet.h"

namespace qmcplusplus {

  /** constructor for EGOSet
   * @param norb number of orbitals for the EGOSet
   * @param k list of unique k points in Cartesian coordinate excluding gamma
   * @param k2 k2[i]=dot(k[i],k[i])
   */
  RealEGOSet::RealEGOSet(const vector<PosType>& k, const vector<RealType>& k2): K(k),mK2(k2) {
    KptMax=k.size();
  }

  ElectronGasOrbitalBuilder::ElectronGasOrbitalBuilder(ParticleSet& els, TrialWaveFunction& psi):
    OrbitalBuilderBase(els,psi)
  { 
  }   


  bool ElectronGasOrbitalBuilder::put(xmlNodePtr cur){

    int nc=0;
    PosType twist(0.0);
    OhmmsAttributeSet aAttrib;
    aAttrib.add(nc,"shell");
    aAttrib.add(twist,"twist");
    aAttrib.put(cur);

    typedef DiracDeterminant<RealEGOSet>  Det_t;
    typedef SlaterDeterminant<RealEGOSet> SlaterDeterminant_t;

    int nat=targetPtcl.getTotalNum();
    int nup=nat/2;

    HEGGrid<RealType,OHMMS_DIM> egGrid(targetPtcl.Lattice);

    if(nc == 0) nc = egGrid.getShellIndex(nup);

    if(nc<0)
    {
      app_error() << "  HEG Invalid Shell." << endl;
      APP_ABORT("ElectronGasOrbitalBuilder::put");
    }

    if(nup!=egGrid.getNumberOfKpoints(nc)) 
    {
      app_error() << "  The number of particles does not match to the shell." << endl;
      app_error() << "  Suggested values for the number of particles " << endl;
      app_error() << "   " << 2*egGrid.getNumberOfKpoints(nc) << " for shell "<< nc << endl;
      app_error() << "   " << 2*egGrid.getNumberOfKpoints(nc-1) << " for shell "<< nc-1 << endl;
      APP_ABORT("ElectronGasOrbitalBuilder::put");
      return false;
    }

    int nkpts=(nup-1)/2;

    //create a E(lectron)G(as)O(rbital)Set
    egGrid.createGrid(nc,nkpts);
    RealEGOSet* psi=new RealEGOSet(egGrid.kpt,egGrid.mk2); 

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
    targetPsi.addOrbital(sdet,"SlaterDet");

    return true;
  }
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
