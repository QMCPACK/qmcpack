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

#include "QMCWaveFunctions/LocalizedBasisSet.h"
#include "QMCWaveFunctions/MolecularOrbitals/AtomicBasisBuilder.h"
#include "QMCWaveFunctions/LCOrbitalSet.h"
#include "Utilities/ProgressReportEngine.h"

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
    typedef LocalizedBasisSet<COT>   ThisBasisSetType;


    /** constructor
     * \param els reference to the electrons
     * \param ions reference to the ions
     */
    MolecularBasisBuilder(ParticleSet& els, ParticleSet& ions):
      targetPtcl(els), sourcePtcl(ions), thisBasisSet(0) 
      { 
        ClassName="MolecularBasisBuilder";
      }   

    bool put(xmlNodePtr cur) 
    {

      if(myBasisSet) return true;

      ReportEngine PRE(ClassName,"put(xmlNodePtr)");

      PRE.echo(cur);

      //create the BasisSetType
      thisBasisSet = new ThisBasisSetType(sourcePtcl,targetPtcl);

      //create the basis set
      //go thru the tree
      cur = cur->xmlChildrenNode;
      while(cur!=NULL) {
        string cname((const char*)(cur->name));
        if(cname == "atomicBasisSet") {
          const xmlChar* eptr=xmlGetProp(cur,(const xmlChar*)"elementType");
          if(eptr == NULL) 
          {
            //fatal error
            PRE.error("Missing elementType attribute of atomicBasisSet.",true);
          }
          string elementType((const char*)eptr);

          map<string,BasisSetBuilder*>::iterator it = aoBuilders.find(elementType); 
          if(it == aoBuilders.end()) {
            AtomicBasisBuilder<RFB>* any = new AtomicBasisBuilder<RFB>(elementType);
            any->setReportLevel(ReportLevel);
            any->put(cur);
            COT* aoBasis= any->createAOSet(cur);
            if(aoBasis) { //add the new atomic basis to the basis set
              int activeCenter =sourcePtcl.getSpeciesSet().findSpecies(elementType);
              thisBasisSet->add(activeCenter, aoBasis);
            }
            aoBuilders[elementType]=any;
          } 
          else 
          {
            PRE.warning("Species "+elementType+" is already initialized. Ignore the input.");
          }
        }
        cur = cur->next;
      }

      //resize the basis set
      thisBasisSet->setBasisSetSize(-1);
      myBasisSet=thisBasisSet;

      return true;
    }

    SPOSetBase* createSPOSet(xmlNodePtr cur) {

      ReportEngine PRE(ClassName,"createSPO(xmlNodePtr)");

      SPOSetBase *lcos=0;
      cur = cur->xmlChildrenNode;
      while(cur!=NULL) {
        string cname((const char*)(cur->name));
        if(cname.find("coeff") < cname.size()) 
        {
          app_log() << "Creating LCOrbitalSet with the input coefficients" << endl;
          lcos= new LCOrbitalSet<ThisBasisSetType,false>(thisBasisSet,ReportLevel);
        }
        cur=cur->next;
      }

      if(lcos==0)
      {
        app_log() << "Creating LCOrbitalSet with the Identity coefficient" << endl;
        lcos = new LCOrbitalSet<ThisBasisSetType,true>(thisBasisSet,ReportLevel);
      }

      return lcos;
    }

  private:
    ///target ParticleSet
    ParticleSet& targetPtcl;
    ///source ParticleSet
    ParticleSet& sourcePtcl;
    ///BasisSet
    ThisBasisSetType* thisBasisSet;
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
