//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_MOBASISBUILDER_H
#define QMCPLUSPLUS_MOBASISBUILDER_H

#include "QMCWaveFunctions/LocalizedBasisSet.h"
#include "QMCWaveFunctions/MolecularOrbitals/AtomicBasisBuilder.h"
#include "QMCWaveFunctions/LCOrbitalSet.h"
#include "QMCWaveFunctions/LCOrbitalSetOpt.h"
#include "Utilities/ProgressReportEngine.h"
#include "OhmmsData/AttributeSet.h"
#include "io/hdf_archive.h"
#include "Message/CommOperators.h"
#if QMC_BUILD_LEVEL>2
#include "QMCWaveFunctions/Experimental/LCOrbitalSetWithCorrection.h"
#endif
namespace qmcplusplus
{

/** derived class from BasisSetBuilder
 *
 * Create a basis set of molecular orbital types as defined in MolecularOrbitalBasis
 * with radial wave functions on the radial grids.
 */
template<class RFB>
class MolecularBasisBuilder: public BasisSetBuilder
{

public:

  typedef typename RFB::CenteredOrbitalType COT;
  typedef LocalizedBasisSet<COT>   ThisBasisSetType;


  /** constructor
   * \param els reference to the electrons
   * \param ions reference to the ions
   */
  MolecularBasisBuilder(ParticleSet& els, ParticleSet& ions, bool cusp=false, std::string cusp_info="",std::string MOH5Ref=""):
    targetPtcl(els), sourcePtcl(ions), thisBasisSet(0), cuspCorr(cusp), cuspInfo(cusp_info), h5_path(MOH5Ref)
  {
    ClassName="MolecularBasisBuilder";
  }

  inline bool is_same(const xmlChar* a, const char* b)
  {
    return !strcmp((const char*)a,b);
  }

  bool put(xmlNodePtr cur)
  {
    if(thisBasisSet)
      return true;


    //Reading from XML
    if(h5_path=="")
    {
      ReportEngine PRE(ClassName,"put(xmlNodePtr)");
      PRE.echo(cur);
      //create the BasisSetType
      thisBasisSet = new ThisBasisSetType(sourcePtcl,targetPtcl);
      if(!is_same(cur->name,"basisset"))
      {//heck to handle things like <sposet_builder>
        xmlNodePtr cur1= cur->xmlChildrenNode;
        while(cur1!=NULL)
        {
          if(is_same(cur1->name,"basisset")) cur=cur1;
          cur1=cur1->next;
        }
      }

      //create the basis set
      //go thru the tree
      cur = cur->xmlChildrenNode;
      while(cur!=NULL)
      {
        std::string cname((const char*)(cur->name));
        if(cname == "atomicBasisSet")
        {
          std::string elementType;
          OhmmsAttributeSet att;
          att.add(elementType,"elementType");
          att.put(cur);
          if(elementType.empty())
            PRE.error("Missing elementType attribute of atomicBasisSet.",true);
          std::map<std::string,BasisSetBuilder*>::iterator it = aoBuilders.find(elementType);
          if(it == aoBuilders.end())
          {
            AtomicBasisBuilder<RFB>* any = new AtomicBasisBuilder<RFB>(elementType);
            any->setReportLevel(ReportLevel);
            any->put(cur);
            COT* aoBasis= any->createAOSet(cur);
            if(aoBasis)
            {
              //add the new atomic basis to the basis set
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
    }
    //Reading from H5
    else
    {
      ReportEngine PRE(ClassName,"put(xmlNodePtr)");
      //create the BasisSetType
      thisBasisSet = new ThisBasisSetType(sourcePtcl,targetPtcl);

      app_log()<<"Reading BasisSet from HDF5 file:"<<h5_path<<std::endl;

      int Nb_Elements(0);
      std::string basiset_name;

      //hdf_archive hin(0);
      hdf_archive hin(myComm);
      if(myComm->rank()==0){
          if(!hin.open(h5_path.c_str(),H5F_ACC_RDONLY))
             PRE.error("Could not open H5 file",true);

          if(!hin.push("basisset"))
             PRE.error("Could not open basisset group in H5; Probably Corrupt H5 file",true);

          hin.read(Nb_Elements,"NbElements");
      }

      myComm->bcast(Nb_Elements);

      if(Nb_Elements<1)
          PRE.error("Missing elementType attribute of atomicBasisSet.",true);



      for (int i=0;i<Nb_Elements;i++)
      {
          std::string elementType,dataset;
          std::stringstream tempElem;
          std::string ElemID0="atomicBasisSet",ElemType;
          tempElem<<ElemID0<<i;
          ElemType=tempElem.str();

          if(myComm->rank()==0){
             if(!hin.push(ElemType.c_str()))
                 PRE.error("Could not open  group Containing atomic Basis set in H5; Probably Corrupt H5 file",true);

             if(!hin.read(basiset_name,"name"))
                 PRE.error("Could not find name of  basisset group in H5; Probably Corrupt H5 file",true);


             if(!hin.read(elementType,"elementType"))
                 PRE.error("Could not read elementType in H5; Probably Corrupt H5 file",true);
          }
          myComm->bcast(basiset_name);
          myComm->bcast(elementType);

          std::map<std::string,BasisSetBuilder*>::iterator it = aoBuilders.find(elementType);
          if(it == aoBuilders.end())
          {
            AtomicBasisBuilder<RFB>* any = new AtomicBasisBuilder<RFB>(elementType);
            any->setReportLevel(ReportLevel);
            any->putH5(hin);
            COT* aoBasis= any->createAOSetH5(hin);
            if(aoBasis)
            {
              //add the new atomic basis to the basis set
              int activeCenter =sourcePtcl.getSpeciesSet().findSpecies(elementType);
              thisBasisSet->add(activeCenter, aoBasis);
            }
            aoBuilders[elementType]=any;
          }

          if(myComm->rank()==0)
             hin.pop();
      }


      if(myComm->rank()==0){
        hin.pop();
        hin.pop();
        hin.close();
      }
    }
    //resize the basis set
    thisBasisSet->setBasisSetSize(-1);
    return true;
  }


  SPOSetBase* createSPOSetFromXML(xmlNodePtr cur)
  {
    ReportEngine PRE(ClassName,"createSPO(xmlNodePtr)");
    std::string spo_name(""), id, cusp_file("");
    std::string use_new_opt_class("no");
    OhmmsAttributeSet spoAttrib;
    spoAttrib.add (spo_name, "name");
    spoAttrib.add (id, "id");
    spoAttrib.add (cusp_file, "cuspInfo");
    spoAttrib.add (use_new_opt_class, "optimize");
    spoAttrib.put(cur);
    SPOSetBase *lcos=0;
    cur = cur->xmlChildrenNode;
    while(cur!=NULL)
    {
      std::string cname((const char*)(cur->name));
      if(cname.find("coeff") < cname.size())
      {
        std::string algorithm("");
        OhmmsAttributeSet coeffAttrib;
        coeffAttrib.add (algorithm, "algorithm");
        coeffAttrib.put(cur);
#if QMC_BUILD_LEVEL>2
        if(cuspCorr)
        {
          app_log() << "Creating LCOrbitalSetWithCorrection with the input coefficients" << std::endl;
          std::string tmp = cuspInfo;
          if(cusp_file != "")
            tmp=cusp_file;
          lcos= new LCOrbitalSetWithCorrection<ThisBasisSetType,false>(thisBasisSet,&targetPtcl,&sourcePtcl,ReportLevel,0.1,tmp,algorithm);
// mmorales:
// this is a small hack to allow the cusp correction to work
// but it should be fixed, all basisset/sposet objects should always be named
          if(spo_name != "")
            lcos->objectName=spo_name;
          else
            lcos->objectName=id;
        }
        else
#endif
        {
          if ( use_new_opt_class == "yes" ) {
            app_log() << "Creating LCOrbitalSetOpt with the input coefficients" << std::endl;
            lcos= new LCOrbitalSetOpt<ThisBasisSetType>(thisBasisSet,ReportLevel);
            if(spo_name != "")
              lcos->objectName=spo_name;
            else
              throw std::runtime_error("LCOrbitalSetOpt spo set must have a name");
          } else {
            app_log() << "Creating LCOrbitalSet with the input coefficients" << std::endl;
            lcos= new LCOrbitalSet<ThisBasisSetType,false>(thisBasisSet,ReportLevel,algorithm);
          }
        }
      }
      cur=cur->next;
    }
    if(lcos==0)
    {
      app_log() << "Creating LCOrbitalSet with the Identity coefficient" << std::endl;
      lcos = new LCOrbitalSet<ThisBasisSetType,true>(thisBasisSet,ReportLevel);
    }
    lcos->myComm=myComm;
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
  std::map<std::string,BasisSetBuilder*> aoBuilders;
  ///apply cusp correction to molecular orbitals
  bool cuspCorr;
  std::string cuspInfo;
  ///read wf info from hdf5
  std::string h5_path;
};
}
#endif
