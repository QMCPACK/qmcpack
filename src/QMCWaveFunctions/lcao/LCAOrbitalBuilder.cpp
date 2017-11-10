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

#include "OhmmsData/AttributeSet.h"
#include <QMCWaveFunctions/SPOSetBase.h>
#include <QMCWaveFunctions/lcao/NGFunctor.h>
#include <QMCWaveFunctions/lcao/MultiQuinticSpline1D.h>
#include "QMCWaveFunctions/lcao/SoaCartesianTensor.h"
#include "QMCWaveFunctions/lcao/SoaSphericalTensor.h"
#include "QMCWaveFunctions/lcao/SoaAtomicBasisSet.h"
#include "QMCWaveFunctions/lcao/SoaLocalizedBasisSet.h"
#include "QMCWaveFunctions/lcao/LCAOrbitalSet.h"
#include "QMCWaveFunctions/lcao/RadialOrbitalSetBuilder.h"
#include "QMCWaveFunctions/lcao/AOBasisBuilder.h"
#include "QMCWaveFunctions/lcao/LCAOrbitalBuilder.h"
#include "QMCWaveFunctions/lcao/MultiFunctorBuilder.h"
#include "Utilities/ProgressReportEngine.h"

namespace qmcplusplus
{
  /** traits for a localized basis set; used by createBasisSet
   *
   * ROT {0=numuerica;, 1=gto; 2=sto} 
   * SH {0=cartesian, 1=spherical}
   * If too confusing, inroduce enumeration.
   */
  template<typename T, int ROT, int SH> struct ao_traits{};

  /** specialization for numerical-cartesian AO */
  template<typename T>
    struct ao_traits<T,0,0>
    {
      typedef MultiQuinticSpline1D<T>                     radial_type;
      typedef SoaCartesianTensor<T>                       angular_type;
      typedef SoaAtomicBasisSet<radial_type,angular_type> ao_type;
      typedef SoaLocalizedBasisSet<ao_type>               basis_type;

    };

  /** specialization for numerical-spherical AO */
  template<typename T>
    struct ao_traits<T,0,1>
    {
      typedef MultiQuinticSpline1D<T>                     radial_type;
      typedef SoaSphericalTensor<T>                       angular_type;
      typedef SoaAtomicBasisSet<radial_type,angular_type> ao_type;
      typedef SoaLocalizedBasisSet<ao_type>               basis_type;
    };

  /** specialization for GTO-cartesian AO */
  template<typename T>
    struct ao_traits<T,1,0>
    {
      typedef MultiFunctorAdapter<GaussianCombo<T> >           radial_type;
      typedef SoaCartesianTensor<T>                       angular_type;
      typedef SoaAtomicBasisSet<radial_type,angular_type> ao_type;
      typedef SoaLocalizedBasisSet<ao_type>               basis_type;
    };

  /** specialization for GTO-cartesian AO */
  template<typename T>
    struct ao_traits<T,1,1>
    {
      typedef MultiFunctorAdapter<GaussianCombo<T> >           radial_type;
      typedef SoaSphericalTensor<T>                       angular_type;
      typedef SoaAtomicBasisSet<radial_type,angular_type> ao_type;
      typedef SoaLocalizedBasisSet<ao_type>               basis_type;
    };

  /** specialization for STO-spherical AO */
  template<typename T>
    struct ao_traits<T,2,1>
    {
      typedef MultiFunctorAdapter<SlaterCombo<T> >             radial_type;
      typedef SoaSphericalTensor<T>                       angular_type;
      typedef SoaAtomicBasisSet<radial_type,angular_type> ao_type;
      typedef SoaLocalizedBasisSet<ao_type>               basis_type;
    };


  inline bool is_same(const xmlChar* a, const char* b)
  {
    return !strcmp((const char*)a,b);
  }


  LCAOrbitalBuilder::LCAOrbitalBuilder(ParticleSet& els, ParticleSet& ions, xmlNodePtr cur) 
    : targetPtcl(els), sourcePtcl(ions), myBasisSet(nullptr)
  {
    ClassName="LCAOrbitalBuilder";
    ReportEngine PRE(ClassName,"createBasisSet");

    std::string keyOpt("NMO"); // Numerical Molecular Orbital
    std::string transformOpt("yes"); // Numerical Molecular Orbital
    std::string cuspC("no");  // cusp correction
    cuspInfo="";
    //std::string cuspInfo("");  // file with precalculated cusp correction info
    OhmmsAttributeSet aAttrib;
    aAttrib.add(keyOpt,"keyword");
    aAttrib.add(keyOpt,"key");
    aAttrib.add(transformOpt,"transform");
    aAttrib.add(cuspC,"cuspCorrection");
    aAttrib.add(cuspInfo,"cuspInfo");
    if(cur != NULL) aAttrib.put(cur);
    
    radialOrbType=-1;
    if (transformOpt == "yes")
      radialOrbType=0;
    else
    {
      if(keyOpt=="GTO") radialOrbType=1;
      if(keyOpt=="STO") radialOrbType=2;
    }

    if(radialOrbType<0)
      PRE.error("Unknown radial function for LCAO orbitals. Specify keyword=\"NMO/GTO/STO\" .",true);
  }

  LCAOrbitalBuilder::~LCAOrbitalBuilder()
  {
    //properly cleanup
  }

  bool LCAOrbitalBuilder::put(xmlNodePtr cur)
  {
    if(myBasisSet != nullptr) return true;

    ReportEngine PRE(ClassName,"put(xmlNodePtr)");

    if(!is_same(cur->name,"basisset"))
    {//heck to handle things like <sposet_builder>
      xmlNodePtr cur1= cur->xmlChildrenNode;
      while(cur1!=NULL)
      {
        if(is_same(cur1->name,"basisset")) cur=cur1;
        cur1=cur1->next;
      }
    }

    int ylm=-1;
    {
      xmlNodePtr cur1= cur->xmlChildrenNode;
      while(cur1!=NULL && ylm<0)
      {
        if(is_same(cur1->name,"atomicBasisSet"))
        {
          std::string sph;
          OhmmsAttributeSet att;
          att.add(sph,"angular");
          att.put(cur1);
          ylm=(sph=="cartesian")?0:1;
        }
        cur1=cur1->next;
      }
    }

    if(ylm<0)
      PRE.error("Missing angular attribute of atomicBasisSet.",true);


    /** process atomicBasisSet per ion species */
    switch(radialOrbType)
    {
      case(0): //numerical
        app_log() << "  LCAO: SoaAtomicBasisSet<MultiQuintic,"<<ylm<<">" << std::endl;;
        if(ylm) 
          myBasisSet=createBasisSet<0,1>(cur);
        else
          myBasisSet=createBasisSet<0,0>(cur);
        break;
      case(1): //gto
        app_log() << "  LCAO: SoaAtomicBasisSet<MultiGTO,"<<ylm<<">" << std::endl;;
        if(ylm) 
          myBasisSet=createBasisSet<1,1>(cur);
        else
          myBasisSet=createBasisSet<1,0>(cur);
        break;
      case(2): //sto
        app_log() << "  LCAO: SoaAtomicBasisSet<MultiSTO,"<<ylm<<">" << std::endl;;
        myBasisSet=createBasisSet<2,1>(cur);
        break;
      default:
        PRE.error("Cannot construct SoaAtomicBasisSet<ROT,YLM>.",true);
        break;
    }
    return true;
  }

  template<int I, int J>
   LCAOrbitalBuilder::BasisSet_t*
   LCAOrbitalBuilder::createBasisSet(xmlNodePtr cur)
  {

    ReportEngine PRE(ClassName,"createBasisSet(xmlNodePtr)");

    typedef typename ao_traits<RealType,I,J>::ao_type    ao_type;
    typedef typename ao_traits<RealType,I,J>::basis_type basis_type;

    basis_type* mBasisSet=new basis_type(sourcePtcl,targetPtcl);

    //keep the builder local
    std::map<std::string,BasisSetBuilder*> aoBuilders;

    /** process atomicBasisSet per ion species */
    cur = cur->xmlChildrenNode;
    while(cur!=NULL) //loop over unique ioons
    {
      std::string cname((const char*)(cur->name));

      if(cname == "atomicBasisSet")
      {
        std::string elementType;
        std::string sph;
        OhmmsAttributeSet att;
        att.add(elementType,"elementType");
        att.put(cur);

        if(elementType.empty())
          PRE.error("Missing elementType attribute of atomicBasisSet.",true);

        std::map<std::string,BasisSetBuilder*>::iterator it = aoBuilders.find(elementType);
        if(it == aoBuilders.end())
        {
          AOBasisBuilder<ao_type>* any = new AOBasisBuilder<ao_type>(elementType);
          any->setReportLevel(ReportLevel);
          any->put(cur);
          ao_type* aoBasis= any->createAOSet(cur);
          if(aoBasis)
          {
            //add the new atomic basis to the basis set
            int activeCenter =sourcePtcl.getSpeciesSet().findSpecies(elementType);
            mBasisSet->add(activeCenter, aoBasis);
          }
          aoBuilders[elementType]=any;
        }
      }
      cur = cur->next;
    } // done with basis set

    { //cleanup basisset builder
      std::map<std::string,BasisSetBuilder*>::iterator itX=aoBuilders.begin();
      while(itX!=aoBuilders.end())
      {
        delete (*itX).second;
        ++itX;
      }
    }

    mBasisSet->setBasisSetSize(-1);
    return mBasisSet;
  }

  SPOSetBase* LCAOrbitalBuilder::createSPOSetFromXML(xmlNodePtr cur)
  {
    ReportEngine PRE(ClassName,"createSPO(xmlNodePtr)");
    std::string spo_name(""), id, cusp_file("");
    OhmmsAttributeSet spoAttrib;
    spoAttrib.add (spo_name, "name");
    spoAttrib.add (id, "id");
    spoAttrib.add (cusp_file, "cuspInfo");
    spoAttrib.put(cur);

    SPOSetBase *lcos=new LCAOrbitalSet(myBasisSet,ReportLevel);

    //@TODO: add cusp condition
    
    return lcos;
  }
}

