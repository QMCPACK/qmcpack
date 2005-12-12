//////////////////////////////////////////////////////////////////
// (c) Copyright 2003  by Jeongnim Kim
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
#include "QMCWaveFunctions/Fermion/SlaterDetBuilder.h"
#include "QMCWaveFunctions/Fermion/DiracDeterminantT.h"
#include "QMCWaveFunctions/MolecularOrbitals/GridMolecularOrbitals.h"
#include "QMCWaveFunctions/LCOrbitals.h"
#include "OhmmsData/AttributeSet.h"

namespace qmcplusplus {

  SlaterDetBuilder::SlaterDetBuilder(ParticleSet& p, 
      TrialWaveFunction& psi, PtclPoolType& psets): 
      OrbitalBuilderBase(p,psi), ptclPool(psets),spoBuilder(0){
    } 

  bool SlaterDetBuilder::put(xmlNodePtr cur) {
    ///vector of Slater determinants 
    std::vector<SlaterDeterminant_t*> slaterdets;
    ///vector of coefficients of the Slater determinants
    std::vector<RealType> sdet_coeff;
    int nvar(targetPsi.VarList.size());
    int is=0, first=0;
    cur = cur->xmlChildrenNode;
    while(cur != NULL) {
      string cname((const char*)(cur->name));
      if(cname == basisset_tag) {
        //create the basis set
        createBasisSet(cur);
      } else if(cname == sd_tag) {
        first = 0;
        //add a new SlaterDeterminant
        slaterdets.push_back(new SlaterDeterminant_t);
        sdet_coeff.push_back(1.0);

        xmlNodePtr tcur = cur->xmlChildrenNode;
        while(tcur != NULL) {
          string tname((const char*)(tcur->name));
          if(tname == param_tag) {
            putContent(sdet_coeff[is],tcur);
          } else if(tname == det_tag) {
            //add the DiracDeterminant to the SlaterDeterminant
            Det_t *adet = createDeterminant(tcur,is,first);
            slaterdets[is]->add(adet);
            first += adet->cols();
          }
          tcur = tcur->next;
        }
        is++;
      }
      cur = cur->next;
    }

    bool optimizeit=(targetPsi.VarList.size()>nvar);
    if(optimizeit) {
      WARNMSG("Slater determinants will be optimized")
    }

    BasisSet->resize(targetPtcl.getTotalNum());

    //distable multi-determinant for now
 //   if(slaterdets.size() > 1) {
 //     MultiSlaterDet *multidet= new MultiSlaterDet;
 //     for(int i=0; i<slaterdets.size(); i++) {
 //       slaterdets[i]->setOptimizable(optimizeit);
 //       multidet->add(slaterdets[i],sdet_coeff[i]);
 //     }
 //     multidet->setOptimizable(optimizeit);
 //     //add a MultiDeterminant to the trial wavefuntion
 //     targetPsi.add(multidet);
 //   } else {
      targetPsi.add(slaterdets[0]);
      if(targetPsi.VarList.size()>nvar) slaterdets[0]->setOptimizable(true);
 //   }
    return true;
  }

  void SlaterDetBuilder::createBasisSet(xmlNodePtr cur) {

    if(spoBuilder ==0) {
      basisName="mo";

      string source("i");
      ParticleSet* ions=0;
      const xmlChar*  a=xmlGetProp(cur,(const xmlChar*)"source");
      if(a) { source = (const char*)a;}
      PtclPoolType::iterator pit = ptclPool.find(source);
      if(pit != ptclPool.end()) ions = (*pit).second;

      GridMolecularOrbitals *gmoBuilder = new GridMolecularOrbitals(targetPtcl,targetPsi,*ions);
      typedef GridMolecularOrbitals::BasisSetType BasisSetType;
      BasisSet = new BasisSetProxy<BasisSetType>(gmoBuilder->addBasisSet(cur));
      spoBuilder=gmoBuilder;
      SPOSetID=0;
    }
  }

  SlaterDetBuilder::Det_t*
  SlaterDetBuilder::createDeterminant(xmlNodePtr cur, int is, int first) {

    string abasis(basisName);
    string detname("invalid"), refname("invalid");
    OhmmsAttributeSet aAttrib;
    aAttrib.add(abasis,basisset_tag);
    aAttrib.add(detname,"id");
    aAttrib.add(refname,"ref");
    aAttrib.put(cur);

    typedef GridMolecularOrbitals::BasisSetType BasisSetType;
    typedef BasisSetProxy<BasisSetType> BasisSetProxyType;
    typedef LCOrbitals<BasisSetType>            SPOSetType;

    SPOSetType *psi=0;
    if(refname == "invalid") { //create one and use detname
      if(detname =="invalid") { //no id is given, assign one
        detname="det";
        ostringstream idassigned(detname);
        idassigned << is;
      }
      BasisSetProxyType* bs = dynamic_cast<BasisSetProxyType*>(BasisSet);
      psi = new SPOSetType(bs->basisRef,is);
      psi->put(cur);
      psi->setName(detname);
      targetPsi.addSPOSet(psi);
    } else {
      if(targetPsi.hasSPOSet(refname)) {
        psi = dynamic_cast<SPOSetType*>(targetPsi.getSPOSet(refname));
      }
    }

    Det_t *adet = new DiracDeterminantT<SPOSetType>(*psi,first);
    adet->set(first,psi->numOrbitals());

    return adet;
  }
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
