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
#include "Utilities/ProgressReportEngine.h"
#include "Numerics/HDFSTLAttrib.h"
#include "OhmmsData/HDFStringAttrib.h"
#include "OhmmsData/AttributeSet.h"
#include "io/hdf_archive.h"
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
  MolecularBasisBuilder(ParticleSet& els, ParticleSet& ions, bool cusp=false, std::string cusp_info="",bool H5Ref=false,std::string MOH5Ref=""):
    targetPtcl(els), sourcePtcl(ions), thisBasisSet(0),cuspCorr(cusp),cuspInfo(cusp_info),h5Ref(H5Ref),h5_path(MOH5Ref)
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
    if(h5Ref==false)
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
      
      std::string basiset_name;
      hid_t h_file,basisset,atomicBasisSet;


      h_file= H5Fopen(h5_path.c_str(),H5F_ACC_RDWR,H5P_DEFAULT);
      if(h_file<0){
         std::cerr<<"Could not open H5 file"<<std::endl;
         exit(0);
      }
      basisset= H5Gopen(h_file,"basisset");
      if(basisset<0){
         std::cerr<<"Could not open basisset group in H5; Probably Corrupt H5 file"<<std::endl;
         exit(0);
      }
      atomicBasisSet= H5Gopen(basisset,"atomicBasisSet");
      if(atomicBasisSet<0){
         std::cerr<<"Could not open atomicBasisSet group in H5; Probably Corrupt H5 file"<<std::endl;
         exit(0);
      }
 
      HDFAttribIO<std::string> basisname(basiset_name); 
      basisname.read(basisset,"name");

      int Nb_Elements(0);
      HDFAttribIO<int> NumElement(Nb_Elements); 
      NumElement.read(atomicBasisSet,"NumElemets");
  
      if(Nb_Elements<1)
          PRE.error("Missing elementType attribute of atomicBasisSet.",true);


      
      for (int i=0;i<Nb_Elements;i++)
      {
          hid_t EleTycBasisSet;  
          std::string elementType,dataset;
          

          std::stringstream tempElem;                                              
          std::string ElemID0="Element",ElemType;
          tempElem<<ElemID0<<i;
          ElemType=tempElem.str();

          EleTycBasisSet= H5Gopen(atomicBasisSet,ElemType.c_str());
          if(atomicBasisSet<0){
            std::cerr<<"Could not open  group EleTycBasisSet in H5; Probably Corrupt H5 file"<<std::endl;
            exit(0);
          }

          HDFAttribIO<std::string> Elementname(elementType);
          Elementname.read(EleTycBasisSet,"elementType");
           
          app_log()<<"BasisSet name: "<<basiset_name<<"    Number of Atoms="<<Nb_Elements<<"     elementType="<<elementType<<std::endl;

          std::map<std::string,BasisSetBuilder*>::iterator it = aoBuilders.find(elementType);
          if(it == aoBuilders.end())
          {
            AtomicBasisBuilder<RFB>* any = new AtomicBasisBuilder<RFB>(elementType);
            any->setReportLevel(ReportLevel);
            any->putH5(EleTycBasisSet);
            COT* aoBasis= any->createAOSetH5(EleTycBasisSet);
            if(aoBasis)
            {
              //add the new atomic basis to the basis set
              int activeCenter =sourcePtcl.getSpeciesSet().findSpecies(elementType);
              thisBasisSet->add(activeCenter, aoBasis);
            }
            aoBuilders[elementType]=any;
          }
          H5Gclose(EleTycBasisSet); 
      }

      H5Gclose(atomicBasisSet);
      H5Gclose(basisset);
      H5Fclose(h_file);
    }
    //resize the basis set
    thisBasisSet->setBasisSetSize(-1);
    return true;
  }


  SPOSetBase* createSPOSetFromXML(xmlNodePtr cur)
  {
    ReportEngine PRE(ClassName,"createSPO(xmlNodePtr)");
    std::string spo_name(""), id, cusp_file("");
    OhmmsAttributeSet spoAttrib;
    spoAttrib.add (spo_name, "name");
    spoAttrib.add (id, "id");
    spoAttrib.add (cusp_file, "cuspInfo");
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
          app_log() << "Creating LCOrbitalSet with the input coefficients" << std::endl;
          lcos= new LCOrbitalSet<ThisBasisSetType,false>(thisBasisSet,ReportLevel,algorithm);
        }
      }
      cur=cur->next;
    }
    if(lcos==0)
    {
      app_log() << "Creating LCOrbitalSet with the Identity coefficient" << std::endl;
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
  std::map<std::string,BasisSetBuilder*> aoBuilders;
  ///apply cusp correction to molecular orbitals
  bool cuspCorr;
  std::string cuspInfo;
  ///read wf info from hdf5
  bool h5Ref;
  std::string h5_path;
};
}
#endif
