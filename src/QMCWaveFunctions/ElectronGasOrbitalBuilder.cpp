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
#include "QMCWaveFunctions/ElectronGasOrbitalBuilder.h"
#include "QMCWaveFunctions/SlaterDeterminant.h"

namespace qmcplusplus {

  /** constructor for EGOSet
   * @param norb number of orbitals for the EGOSet
   * @param k list of unique k points in Cartesian coordinate excluding gamma
   * @param k2 k2[i]=dot(k[i],k[i])
   */
  ElectronGasOrbitalBuilder::EGOSet::EGOSet(const vector<PosType>& k, const vector<RealType>& k2): K(k),mK2(k2) {
    KptMax=k.size();
  }

  ElectronGasOrbitalBuilder::ElectronGasOrbitalBuilder(ParticleSet& els, TrialWaveFunction& psi):
    OrbitalBuilderBase(els,psi)
  { 
  }   


  bool ElectronGasOrbitalBuilder::put(xmlNodePtr cur){

    int nc=0;
    const xmlChar* nc_ptr=xmlGetProp(cur,(const xmlChar*)"shell");
    if(nc_ptr) {
      nc = atoi((const char*)nc_ptr);
    }

    typedef DiracDeterminant<EGOSet>  Det_t;
    typedef SlaterDeterminant<EGOSet> SlaterDeterminant_t;

    int nat=targetPtcl.getTotalNum();
    int nup=nat/2;
    if(nc == 0) {
      nc = static_cast<int>(pow(static_cast<double>(nup),1.0/3.0))/2+1;
    }

    cout << "   The number of shells " << nc << endl;
    map<int,vector<PosType>*> rs;
    int first_ix2, first_ix3; 
    for(int ix1=0; ix1<=nc; ix1++) {
      if(ix1 == 0) first_ix2=0;
      else         first_ix2=-nc;
      for(int ix2=first_ix2; ix2<=nc; ix2++) {
        if(ix1 == 0 && ix2 == 0) first_ix3=1;
        else                     first_ix3=-nc;
        for(int ix3=first_ix3; ix3<=nc; ix3++) {
          int ih=ix1*ix1+ix2*ix2+ix3*ix3;
          std::map<int,vector<PosType>*>::iterator it = rs.find(ih);
          if(it == rs.end()) {
            vector<PosType>* ns = new vector<PosType>;
            ns->push_back(PosType(ix1,ix2,ix3));
            rs[ih] = ns;
          } else {
            (*it).second->push_back(PosType(ix1,ix2,ix3));
          }
        }
      }
    }

    int nkpts=(nup-1)/2;

    vector<PosType> kpt(nkpts);
    vector<RealType> mk2(nkpts);

    int ikpt=0;
    map<int, vector<PosType>*>::iterator rs_it(rs.begin()), rs_end(rs.end());
    while(ikpt<nkpts && rs_it != rs_end) {
      vector<PosType>::iterator ns_it((*rs_it).second->begin()), ns_end((*rs_it).second->end());
      RealType minus_ksq=-targetPtcl.Lattice.ksq(*ns_it);
      while(ikpt<nkpts && ns_it!=ns_end) {
        kpt[ikpt]=targetPtcl.Lattice.k_cart(*ns_it);
        mk2[ikpt]=minus_ksq;
        ++ikpt;
        ++ns_it;
      }
      ++rs_it;
    }

    cout << "Lattice of electrons " << endl;
    targetPtcl.Lattice.print(cout);

    cout << "Number of kpts " << nkpts << endl;

    cout << "Initial position of the electrons " << endl;
    for(int iat=0; iat<nat; iat++) {
      cout << targetPtcl.R[iat] <<endl;
    }
    cout << "List of kpoints (half-sphere) " << endl;
    for(int ik=0; ik<kpt.size(); ik++) {
      cout << ik << " " << kpt[ik] << " " << mk2[ik] << endl;
    }

    //create a E(lectron)G(as)O(rbital)Set
    EGOSet* psi=new EGOSet(kpt,mk2); 

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
