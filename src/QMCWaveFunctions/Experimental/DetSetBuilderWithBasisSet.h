//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#ifndef QMCPLUSPLUS_DETSET_BUILDER_WITH_BASISSET_H
#define QMCPLUSPLUS_DETSET_BUILDER_WITH_BASISSET_H

#include "Numerics/LibxmlNumericIO.h"
#include "QMCWaveFunctions/OrbitalBuilderBase.h"
#include "QMCWaveFunctions/SlaterDeterminant.h"
#include "QMCWaveFunctions/MultiSlaterDeterminant.h"
#include "QMCWaveFunctions/LCOrbitals.h"

namespace qmcplusplus
{

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
struct DetSetBuilderWithBasisSet: public OrbitalBuilderBase
{

  BasisBuilderT& builder_ref;
  int NumPtcl;
  std::string BasisName;

  /** constructor
   *@param p particleset whose positions defines the wave function
   *@param psi trial wavefuntion to which determinant terms are added
   *@param abuilder a BasisBuilderT object, provides addBasisSet and typedefs
   */
  DetSetBuilderWithBasisSet(ParticleSet& p, TrialWaveFunction& psi, BasisBuilderT& abuilder):
    OrbitalBuilderBase(p,psi), builder_ref(abuilder), BasisName("mo")
  {
  }

  /** process the current xml node to create single-particle orbital
   *@param cur xmlNodePtr to be processed
   *@return true when the size of determinants is positive.
   */
  bool put(xmlNodePtr cur)
  {
    typedef typename BasisBuilderT::BasisSetType BasisSetType;
    typedef LCOrbitals<BasisSetType>             SPOSetType;
    typedef DiracDeterminant<SPOSetType>         Det_t;
    typedef SlaterDeterminant<SPOSetType>        SlaterDeterminant_t;
    ///vector of Slater determinants
    std::vector<SlaterDeterminant_t*> slaterdets;
    std::vector<std::string> sdet_coeff_id;
    std::vector<RealType> sdet_coeff;
    ///vector of coefficients of the Slater determinants
    //std::vector<RealType> sdet_coeff;
    ///pointer to the basis set
    BasisSetType *basisSet =0;
    int nvar(targetPsi.VarList.size());
    int is=0, first=0;
    int detCounter = 0;
    std::string cname, tname;
    cur = cur->xmlChildrenNode;
    while(cur != NULL)
    {
      getNodeName(cname,cur);
      if(cname == basisset_tag)
      {
        //call the BasisSet builder
        basisSet = builder_ref.addBasisSet(cur);
        if(!basisSet)
          return 0;
      }
      else
        if(cname == sd_tag)
        {
          first = 0;
          //add a new SlaterDeterminant
          slaterdets.push_back(new SlaterDeterminant_t);
          std::string multiCoeff_name("detC");
          sdet_coeff.push_back(1.0);
          xmlNodePtr tcur = cur->xmlChildrenNode;
          while(tcur != NULL)
          {
            //string tname((const char*)(tcur->name));
            getNodeName(tname,tcur);
            if(tname == param_tag)
            {
              const xmlChar* cptr=xmlGetProp(tcur,(const xmlChar*)"id");
              if(cptr != NULL)
                multiCoeff_name=(const char*)cptr;
              putContent(sdet_coeff[is],tcur);
            }
            else
              if(tname == det_tag)
              {
                std::string basisName(BasisName);
                std::string detname("invalid"), refname("invalid");
                OhmmsAttributeSet aAttrib;
                aAttrib.add(basisName,basisset_tag);
                aAttrib.add(detname,"id");
                aAttrib.add(refname,"ref");
                aAttrib.put(tcur);
                if(basisName == BasisName)
                  //only when the basisset matches
                {
                  SPOSetType* psi=0;
                  if(refname == "invalid")
                    //create one and use detname
                  {
                    if(detname =="invalid")
                      //no id is given, assign one
                    {
                      detname="det";
                      std::ostringstream idassigned(detname);
                      idassigned << is;
                    }
                    psi = new SPOSetType(basisSet,detCounter);
                    psi->put(tcur);
                    psi->setName(detname);
                    targetPsi.addSPOSet(psi);
                  }
                  else
                  {
                    if(targetPsi.hasSPOSet(refname))
                    {
                      psi = dynamic_cast<SPOSetType*>(targetPsi.getSPOSet(refname));
                    }
                  }
                  Det_t *adet = new Det_t(*psi,first);
                  adet->set(first,psi->numOrbitals());
                  XMLReport("Add a determinant to the SlaterDeterminant for particles: " << first << " -> " << first+psi->numOrbitals())
                  //add the DiracDeterminant to the SlaterDeterminant
                  slaterdets[is]->add(adet);
                  first += psi->numOrbitals();
                  //increment detCounter so that only one SPOSetType object handles the BasisSet
                  detCounter++;
                }
              }
            tcur = tcur->next;
            sdet_coeff_id.push_back(multiCoeff_name);
          }//cname == sd_tag
          is++;
        }
      cur = cur->next;
    }
    bool optimizeit=(targetPsi.VarList.size()>nvar);
    if(optimizeit)
    {
      WARNMSG("Slater determinants will be optimized")
    }
    //set the basis set
    XMLReport("Assigning BasisSet to the first SlaterDeterminant")
    slaterdets[0]->setBasisSet(basisSet);
    if(slaterdets.size() > 1)
    {
      XMLReport("Creating a multi-determinant wavefunction")
      MultiSlaterDeterminant<SPOSetType>
      *multidet= new MultiSlaterDeterminant<SPOSetType>;
      for(int i=0; i<slaterdets.size(); i++)
      {
        slaterdets[i]->setOptimizable(optimizeit);
        multidet->add(slaterdets[i],sdet_coeff[i],sdet_coeff_id[i],targetPsi.VarList);
        app_log() << "   MultiSlaterDeterminant coeff name=" << sdet_coeff_id[i]
                  << " value=" << sdet_coeff[i] << std::endl;
      }
      multidet->setOptimizable(true);
      //add a MultiDeterminant to the trial wavefuntion
      targetPsi.addOrbital(multidet);
    }
    else
    {
      XMLReport("Creating a SlaterDeterminant wavefunction")
      //add a SlaterDeterminant to the trial wavefuntion
      targetPsi.addOrbital(slaterdets[0]);
      if(targetPsi.VarList.size()>nvar)
        slaterdets[0]->setOptimizable(true);
    }
    NumPtcl=first;
    return (NumPtcl> 0);
  }
};
}
#endif
