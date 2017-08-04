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
// File created by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//////////////////////////////////////////////////////////////////////////////////////

#include "Utilities/ProgressReportEngine.h"
#include "OhmmsData/AttributeSet.h"
#include "QMCWaveFunctions/lcao/LCAOrbitalBuilder.h"
#include "QMCWaveFunctions/lcao/LCAOrbitalSet.h"
#include "QMCWaveFunctions/lcao/AOBasisBuilder.h"

namespace qmcplusplus
{
  LCAOrbitalBuilder::LCAOrbitalBuilder(ParticleSet& els, ParticleSet& ions, bool cusp, std::string cusp_info):
    targetPtcl(els), sourcePtcl(ions), xyzBasisSet(nullptr), ylmBasisSet(nullptr),cuspCorr(cusp),cuspInfo(cusp_info)
  {
    ClassName="MolecularBasisBuilder";
  }

  bool LCAOrbitalBuilder::put(xmlNodePtr cur)
  {
    if(xyzBsisSet != nullptr || ylmBasisSet != nullptr) return true;

    ReportEngine PRE(ClassName,"put(xmlNodePtr)");
    PRE.echo(cur);

    bool useCartesian=false;

    if(!is_same(cur->name,"basisset"))
    {//heck to handle things like <sposet_builder>
      xmlNodePtr cur1= cur->xmlChildrenNode;
      while(cur1!=NULL)
      {
        if(is_same(cur1->name,"basisset")) cur=cur1;
        cur1=cur1->next;
      }
    }

    /** process atomicBasisSet per ion species */
    cur = cur->xmlChildrenNode;
    while(cur!=NULL)
    {
      std::string cname((const char*)(cur->name));
      if(cname == "atomicBasisSet")
      {
        std::string elementType;
        std::string sph;
        std::string Morder("gaussian");

        OhmmsAttributeSet att;
        att.add(elementType,"elementType");
        att.add(sph,"angular");
        att.add(Morder,"expandYlm");
        att.put(cur);

        if(elementType.empty())
          PRE.error("Missing elementType attribute of atomicBasisSet.",true);
        std::map<std::string,BasisSetBuilder*>::iterator it = aoBuilders.find(elementType);
        if(it == aoBuilders.end())
        {
          if(sph == "cartesian" || Morder == "Gamess")
          {
            if(xyzBasisSet == nullptr)
              xyzBasisSet = new XYZBasisT(sourcePtcl,targetPtcl);
            AOBasisBuilder<XYZCOT>* any = new AOBasisBuilder<XYZCOT>(elementType);
            any->setReportLevel(ReportLevel);
            any->put(cur);
            XYZCOT* aoBasis= any->createAOSet(cur);
            if(aoBasis)
            {
              //add the new atomic basis to the basis set
              int activeCenter =sourcePtcl.getSpeciesSet().findSpecies(elementType);
              xyzBasisSet->add(activeCenter, aoBasis);
            }
            aoBuilders[elementType]=any;
          }
          else
          {
            if(ylmBasisSet == nullptr)
              ylmBasisSet = new YlmBasisT(sourcePtcl,targetPtcl);
            AOBasisBuilder<YlmCOT>* any = new AOBasisBuilder<YlmCOT>(elementType);
            any->setReportLevel(ReportLevel);
            any->put(cur);
            YlmCOT* aoBasis= any->createAOSet(cur);
            if(aoBasis)
            {
              //add the new atomic basis to the basis set
              int activeCenter =sourcePtcl.getSpeciesSet().findSpecies(elementType);
              YlmBasisSet->add(activeCenter, aoBasis);
            }
            aoBuilders[elementType]=any;
          }
        }
      } 
      else
      {
        PRE.warning("Species "+elementType+" is already initialized. Ignore the input.");
      }
      cur = cur->next;
    }
    return true;
  }

  SPOSetBase* LCAOrbtialBuilder::createSPOSetFromXML(xmlNodePtr cur)
  {
    ReportEngine PRE(ClassName,"createSPO(xmlNodePtr)");
    std::string spo_name(""), id, cusp_file("");
    OhmmsAttributeSet spoAttrib;
    spoAttrib.add (spo_name, "name");
    spoAttrib.add (id, "id");
    spoAttrib.add (cusp_file, "cuspInfo");
    spoAttrib.put(cur);
    SPOSetBase *lcos=nullptr;
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
        app_log() << "Creating LCOrbitalSet with the input coefficients" << std::endl;
        if(xyzBasisSet!=nullptr)
        {
          lcos= new LCAOrtbialSet<XYZBasisT>(xyzBasisSet,ReportLevel);
          //take care of the cusp condition
        }
        if(ylmBasisSet!=nullptr)
        {
          lcos= new LCAOrtbialSet<YlmBasisT>(ylmBasisSet,ReportLevel);
          //take care of the cusp condition
        }
        //#if QMC_BUILD_LEVEL>2
        //            if(cuspCorr)
        //            {
        //              app_log() << "Creating LCOrbitalSetWithCorrection with the input coefficients" << std::endl;
        //              std::string tmp = cuspInfo;
        //              if(cusp_file != "")
        //                tmp=cusp_file;
        //              lcos= new LCOrbitalSetWithCorrection<ThisBasisSetType,false>(thisBasisSet,&targetPtcl,&sourcePtcl,ReportLevel,0.1,tmp,algorithm);
        //              // mmorales:
        //              // this is a small hack to allow the cusp correction to work
        //              // but it should be fixed, all basisset/sposet objects should always be named
        //              if(spo_name != "")
        //                lcos->objectName=spo_name;
        //              else
        //                lcos->objectName=id;
        //            }
        //            else
        //#endif
      }
      cur=cur->next;
    }
    if(lcos==0)
    {//rare case
      app_log() << "Creating LCOrbitalSet with the Identity coefficient" << std::endl;
      if(xyzBasisSet!=nullptr)
        lcos= new LCAOrtbialSet<XYZBasisT>(xyzBasisSet,ReportLevel,true);
      if(ylmBasisSet!=nullptr)
        lcos= new LCAOrtbialSet<YlmBasisT>(ylmBasisSet,ReportLevel,true);
      lcos->setIdentity(true);
    }
    return lcos;
  }
}

