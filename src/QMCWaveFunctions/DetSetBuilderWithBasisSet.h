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
#ifndef OHMMS_QMC_DETSET_BUILDER_WITH_BASISSET_H
#define OHMMS_QMC_DETSET_BUILDER_WITH_BASISSET_H

#include "Numerics/LibxmlNumericIO.h"
#include "QMCWaveFunctions/OrbitalBuilderBase.h"
#include "QMCWaveFunctions/SlaterDeterminant.h"
#include "QMCWaveFunctions/MultiSlaterDeterminant.h"
#include "QMCWaveFunctions/LCOrbitals.h"

namespace ohmmsqmc {

  /**Class to add SlaterDeterminant/MultiSlaterDeterminant to a TrialWaveFunction.
   *
   *The template class BasisBuilderT should provide
   * - BasisSetType: the type of BasisSet
   * - BasisSetType* addBasisSet(xmlNodePtr): returns a pointer to the BasisSet
   *
   * Examples of a BasisBuilderT:
   * - MolecularOrbitals/STOMolecularOrbitals.h
   * - MolecularOrbitals/GridMolecularOrbitals.h
   *
   *@note A DiracDeterminant is a single determinant, a SlaterDeterminant is 
   *the product of DiracDeterminants while a MultiDeterminant is a linear
   *combination of SlaterDeterminants
   */
  template<class BasisBuilderT>
  struct DetSetBuilderWithBasisSet: public OrbitalBuilderBase {
  
    BasisBuilderT& builder_ref;
    int NumPtcl;

    DetSetBuilderWithBasisSet(TrialWaveFunction& awfs, 
			      BasisBuilderT& abuilder): 
      OrbitalBuilderBase(awfs),
      builder_ref(abuilder){
    } 
 
    /** process the current xml node to create single-particle orbital
     *@param wfs_ref a trial wavefuntion to which determinant terms are added
     *@param builder a BasisBuilderT object, provides addBasisSet and typedefs
     *@param cur xmlNodePtr to be processed
     *@return the number of the quantum particles for the determinant terms
     */
    bool put(xmlNodePtr cur) {
    
      typedef typename BasisBuilderT::BasisSetType BasisSetType;
      typedef LCOrbitals<BasisSetType>             SPOSetType;
      typedef DiracDeterminant<SPOSetType>         Det_t;
      typedef SlaterDeterminant<SPOSetType>        SlaterDeterminant_t;
      
      ///vector of Slater determinants 
      std::vector<SlaterDeterminant_t*> slaterdets;
      ///vector of coefficients of the Slater determinants
      std::vector<RealType> sdet_coeff;
      ///pointer to the basis set
      BasisSetType *basisSet =NULL;
      
      //std::map<std::string,SPOSetType*> spoSet;

      int nvar(wfs_ref.VarList.size());
      int is=0, first=0;
      int norbs = 0;
      cur = cur->xmlChildrenNode;
      while(cur != NULL) {
	string cname((const char*)(cur->name));
	if(cname == basisset_tag) {
	  //call the BasisSet builder
	  basisSet = builder_ref.addBasisSet(cur);
	  if(!basisSet) return 0;
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
              SPOSetType* psi=0;
              string detname("det");
              const xmlChar* a=xmlGetProp(tcur,(const xmlChar*)"id");
              if(a) {
                detname = (const char*)a;
              } else {
                ostringstream idassigned(detname);
                idassigned << is;
              }

              a=xmlGetProp(tcur,(const xmlChar*)"ref");
              if(a) {
                string detref((const char*)a);
                if(wfs_ref.hasSPOSet(detref)) {
                  psi = dynamic_cast<SPOSetType*>(wfs_ref.getSPOSet(detref));
                }
                //std::map<std::string,SPOSetType*>::iterator dit(spoSet.find(detref));
                //if(dit == spoSet.end()) {
                //  ERRORMSG(detref << " has not be added. Using the first determinant")
                //  psi=(*(spoSet.begin())).second;
                //} else {
                //  XMLReport(detname << " uses " << detref << " determinant.")
                //  psi=(*dit).second;
                //}
              } else {
 	        psi = new SPOSetType(basisSet,norbs);
	        psi->put(tcur);
                psi->setName(detname);
                wfs_ref.addSPOSet(psi);
                //spoSet[detname]=psi;
              }
	      Det_t *adet = new Det_t(*psi,first);
	      adet->set(first,psi->numOrbitals());
	      XMLReport("Add a determinant to the SlaterDeterminant for particles: " 
			<< first << " -> " << first+psi->numOrbitals())
		//add the DiracDeterminant to the SlaterDeterminant
	      slaterdets[is]->add(adet);
	      first += psi->numOrbitals();
	      norbs++;
	    }
	    tcur = tcur->next;
	  }
	  is++;
	}
	cur = cur->next;
      }
      
      bool optimizeit=(wfs_ref.VarList.size()>nvar);
      if(optimizeit) {
        WARNMSG("Slater determinants will be optimized")
      }
      if(slaterdets.size() > 1) {
	XMLReport("Creating a multi-determinant wavefunction")
	MultiSlaterDeterminant<SPOSetType>
	  *multidet= new MultiSlaterDeterminant<SPOSetType>;
	for(int i=0; i<slaterdets.size(); i++) {
	  XMLReport("Coefficient for a SlaterDeterminant " << sdet_coeff[i])
          slaterdets[i]->setOptimizable(optimizeit);
	  multidet->add(slaterdets[i],sdet_coeff[i]);
	}
        multidet->setOptimizable(optimizeit);
	//add a MultiDeterminant to the trial wavefuntion
	wfs_ref.add(multidet);
      } else {
	XMLReport("Creating a SlaterDeterminant wavefunction")
	//add a SlaterDeterminant to the trial wavefuntion
	wfs_ref.add(slaterdets[0]);
        if(wfs_ref.VarList.size()>nvar) slaterdets[0]->setOptimizable(true);
	//add a MultiDeterminant to the trial wavefuntion
      }
      NumPtcl=first; 
      return (NumPtcl> 0);
    }
  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
