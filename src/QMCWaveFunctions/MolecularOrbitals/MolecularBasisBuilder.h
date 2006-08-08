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
#ifndef QMCPLUSPLUS_MOBASISBUILDER_H
#define QMCPLUSPLUS_MOBASISBUILDER_H

#include "QMCWaveFunctions/SphericalOrbitalSet.h"
#include "QMCWaveFunctions/LocalizedBasisSet.h"
#include "QMCWaveFunctions/MolecularOrbitals/AtomicBasisBuilder.h"
#include "QMCWaveFunctions/LCOrbitalSet.h"
#include "Message/Communicate.h"

namespace qmcplusplus {


 /** derived class from BasisSetBuilder 
  *
  * Create a basis set of molecular orbital types as defined in MolecularOrbitalBasis
  * with radial wave functions on the radial grids.
  */
  template<class RFB>
  class MolecularBasisBuilder: public BasisSetBuilder {

  public:

    typedef typename RFB::CenteredOrbitalType COT;
    typedef LocalizedBasisSet<COT>   BasisSetType;

    /** constructor
     * \param els reference to the electrons
     * \param ions reference to the ions
     */
    MolecularBasisBuilder(ParticleSet& els, ParticleSet& ions):
      targetPtcl(els), sourcePtcl(ions), myBasisSet(0) { }   

    bool put(xmlNodePtr cur) {

      if(myBasisSet) return true;

      DistanceTableData* d_table=DistanceTable::add(sourcePtcl,targetPtcl);

      //create the BasisSetType
      myBasisSet = new BasisSetType(sourcePtcl,targetPtcl);

      //create the basis set
      //go thru the tree
      cur = cur->xmlChildrenNode;
      while(cur!=NULL) {
        string cname((const char*)(cur->name));
        if(cname == "atomicBasisSet") {
          const xmlChar* eptr=xmlGetProp(cur,(const xmlChar*)"elementType");
          if(eptr == NULL) {
            app_error() << "   Missing elementType attribute of atomicBasisSet.\n"
              << " Abort at MOBasisBuilder::put " << endl;
            OHMMS::Controller->abort();
          }
          string elementType((const char*)eptr);

          map<string,BasisSetBuilder*>::iterator it = aoBuilders.find(elementType); 
          if(it == aoBuilders.end()) {
            AtomicBasisBuilder<RFB>* any = new AtomicBasisBuilder<RFB>(elementType);
            any->put(cur);
            COT* aoBasis= any->createAOSet(cur);
            if(aoBasis) { //add the new atomic basis to the basis set
              int activeCenter =sourcePtcl.getSpeciesSet().findSpecies(elementType);
              myBasisSet->add(activeCenter, aoBasis);
            }
            aoBuilders[elementType]=any;
          } else {
            WARNMSG("Species " << elementType << " is already initialized. Ignore the input.")
          }
        }
        cur = cur->next;
      }

      //resize the basis set
      myBasisSet->setBasisSetSize(-1);
      return true;
    }

    SPOSetBase* createSPOSet(xmlNodePtr cur) {
      return new LCOrbitalSet<BasisSetType>(myBasisSet);
    }

  private:
    ///target ParticleSet
    ParticleSet& targetPtcl;
    ///source ParticleSet
    ParticleSet& sourcePtcl;
    ///BasisSet
    BasisSetType* myBasisSet;
    ///save AtomiBasisBuilder<RFB>*
    map<string,BasisSetBuilder*> aoBuilders;
  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
