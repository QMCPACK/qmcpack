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
#include "io/hdf_archive.h"
#include "Message/CommOperators.h"
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
    : targetPtcl(els), sourcePtcl(ions), myBasisSet(nullptr), h5_path("")
  {
    ClassName="LCAOrbitalBuilder";
    ReportEngine PRE(ClassName,"createBasisSet");

    std::string keyOpt("NMO"); // Numerical Molecular Orbital
    std::string transformOpt("yes"); // Numerical Molecular Orbital
    std::string cuspC("no");  // cusp correction
    PosType Images(0.0);
    cuspInfo="";  // file with precalculated cusp correction info
    OhmmsAttributeSet aAttrib;
    aAttrib.add(keyOpt,"keyword");
    aAttrib.add(keyOpt,"key");
    aAttrib.add(transformOpt,"transform");
    aAttrib.add(cuspC,"cuspCorrection");
    aAttrib.add(cuspInfo,"cuspInfo");
    aAttrib.add(h5_path,"href");
    aAttrib.add(Images,"PBCimages");
    aAttrib.put(cur);
    PBCImages.resize(3);
    PBCImages[0]=int(Images[0]);
    PBCImages[1]=int(Images[1]);
    PBCImages[2]=int(Images[2]);
     
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
    if(h5_path=="")
      putXML(cur);
    else
      putH5();
    return true;
  }


  bool LCAOrbitalBuilder::putXML(xmlNodePtr cur)
  {
    ReportEngine PRE(ClassName,"putXML(xmlNodePtr)");
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
        app_log() << "  LCAO: SoaAtomicBasisSet<MultiQuintic,"<<ylm<<">" << std::endl;
        if(ylm) 
          myBasisSet=createBasisSet<0,1>(cur);
        else
          myBasisSet=createBasisSet<0,0>(cur);
        break;
      case(1): //gto
        app_log() << "  LCAO: SoaAtomicBasisSet<MultiGTO,"<<ylm<<">" << std::endl;
        if(ylm) 
          myBasisSet=createBasisSet<1,1>(cur);
        else
          myBasisSet=createBasisSet<1,0>(cur);
        break;
      case(2): //sto
        app_log() << "  LCAO: SoaAtomicBasisSet<MultiSTO,"<<ylm<<">" << std::endl;
        myBasisSet=createBasisSet<2,1>(cur);
        break;
      default:
        PRE.error("Cannot construct SoaAtomicBasisSet<ROT,YLM>.",true);
        break;
    }
    return true;
  }

  bool LCAOrbitalBuilder::putH5()
  {
    ReportEngine PRE(ClassName,"putH5()");
    hdf_archive hin(myComm);
    int ylm=-1;
    if(myComm->rank()==0)
    {
      if(!hin.open(h5_path.c_str(),H5F_ACC_RDONLY))
        PRE.error("Could not open H5 file",true);
      if(!hin.push("basisset"))
        PRE.error("Could not open basisset group in H5; Probably Corrupt H5 file",true);

      std::string sph;
      std::string ElemID0="atomicBasisSet0";
      if(!hin.push(ElemID0.c_str()))
        PRE.error("Could not open  group Containing atomic Basis set in H5; Probably Corrupt H5 file",true);
      if(!hin.read(sph,"angular"))
        PRE.error("Could not find name of  basisset group in H5; Probably Corrupt H5 file",true);
      ylm=(sph=="cartesian")?0:1;
      hin.close();
    }

    myComm->bcast(ylm);
    if(ylm<0)
      PRE.error("Missing angular attribute of atomicBasisSet.",true);

    /** process atomicBasisSet per ion species */
    switch(radialOrbType)
    {
      case(0): //numerical
        app_log() << "  LCAO: SoaAtomicBasisSet<MultiQuintic,"<<ylm<<">" << std::endl;
        if(ylm) 
          myBasisSet=createBasisSetH5<0,1>();
        else
          myBasisSet=createBasisSetH5<0,0>();
        break;
      case(1): //gto
        app_log() << "  LCAO: SoaAtomicBasisSet<MultiGTO,"<<ylm<<">" << std::endl;
        if(ylm) 
          myBasisSet=createBasisSetH5<1,1>();
        else
          myBasisSet=createBasisSetH5<1,0>();
        break;
      case(2): //sto
        app_log() << "  LCAO: SoaAtomicBasisSet<MultiSTO,"<<ylm<<">" << std::endl;
        myBasisSet=createBasisSetH5<2,1>();
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
    mBasisSet->setPBCImages(PBCImages);
    return mBasisSet;
  }


  template<int I, int J>
   LCAOrbitalBuilder::BasisSet_t*
   LCAOrbitalBuilder::createBasisSetH5()
  {

    ReportEngine PRE(ClassName,"createBasisSetH5(xmlNodePtr)");

    typedef typename ao_traits<RealType,I,J>::ao_type    ao_type;
    typedef typename ao_traits<RealType,I,J>::basis_type basis_type;

    basis_type* mBasisSet=new basis_type(sourcePtcl,targetPtcl);

    //keep the builder local
    std::map<std::string,BasisSetBuilder*> aoBuilders;
    
    int Nb_Elements(0);
    std::string basiset_name;

    /** process atomicBasisSet per ion species */
    app_log() << "Reading BasisSet from HDF5 file:" << h5_path << std::endl;

    hdf_archive hin(myComm);
    if(myComm->rank()==0)
    {
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

      if(myComm->rank()==0)
      {
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
        AOBasisBuilder<ao_type>* any = new AOBasisBuilder<ao_type>(elementType);
        any->setReportLevel(ReportLevel);
        any->putH5(hin);
        ao_type* aoBasis= any->createAOSetH5(hin);
        if(aoBasis)
        {
          //add the new atomic basis to the basis set
          int activeCenter =sourcePtcl.getSpeciesSet().findSpecies(elementType);
          mBasisSet->add(activeCenter, aoBasis);
        }
        aoBuilders[elementType]=any;
      }

      if(myComm->rank()==0)
        hin.pop();
    }

    if(myComm->rank()==0)
    {
      hin.pop();
      hin.close();
    }

    { //cleanup basisset builder
      std::map<std::string,BasisSetBuilder*>::iterator itX=aoBuilders.begin();
      while(itX!=aoBuilders.end())
      {
        delete (*itX).second;
        ++itX;
      }
    }

    mBasisSet->setBasisSetSize(-1);
    mBasisSet->setPBCImages(PBCImages);
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
    lcos->myComm=myComm;

    //@TODO: add cusp condition
    
    return lcos;
  }


bool SPOSetBase::put(xmlNodePtr cur)
{
  #undef FunctionName
#define FunctionName printf("Calling FunctionName from %s\n",__FUNCTION__);FunctionNameReal
  //Check if HDF5 present
  ReportEngine PRE("SPOSetBase","put(xmlNodePtr)");

  //Special case for sposet hierarchy: go up only once.
  OhmmsAttributeSet locAttrib;
  std::string cur_name;
  locAttrib.add (cur_name, "name");
  locAttrib.put(cur);
  xmlNodePtr curtemp;
  if (cur_name=="spo-up" || cur_name=="spo-dn")
     curtemp=cur->parent;
  else
     curtemp=cur->parent->parent;
     
  std::string MOtype,MOhref;
  bool H5file=false;
  OhmmsAttributeSet H5checkAttrib;
  H5checkAttrib.add(MOtype,"type");
  H5checkAttrib.add(MOhref,"href");
  H5checkAttrib.put(curtemp);
  xmlChar* MOhreftemp;
  if(MOtype=="MolecularOrbital" && MOhref!="")
  {
     MOhreftemp=xmlGetProp(curtemp, (xmlChar*)"href");
     H5file=true;
     PRE.echo(curtemp);
  }

  const char* MOhref2((const char*)MOhreftemp);

  //initialize the number of orbital by the basis set size
  int norb= BasisSetSize;
  std::string debugc("no");
  double orbital_mix_magnitude = 0.0;
  bool PBC=false;
  OhmmsAttributeSet aAttrib;
  aAttrib.add(norb,"orbitals");
  aAttrib.add(norb,"size");
  aAttrib.add(debugc,"debug");
  aAttrib.add(orbital_mix_magnitude, "orbital_mix_magnitude");
  aAttrib.put(cur);
  setOrbitalSetSize(norb);
  xmlNodePtr occ_ptr=NULL;
  xmlNodePtr coeff_ptr=NULL;
  cur = cur->xmlChildrenNode;
  while(cur != NULL)
  {
    std::string cname((const char*)(cur->name));
    if(cname == "occupation")
    {
      occ_ptr=cur;
    }
    else if(cname.find("coeff") < cname.size() || cname == "parameter" || cname == "Var")
    {
      coeff_ptr=cur;
    }
    cur=cur->next;
  }
  if(coeff_ptr == NULL)
  {
    app_log() << "   Using Identity for the LCOrbitalSet " << std::endl;
    return setIdentity(true);
  }
  bool success=putOccupation(occ_ptr);
  if(H5file==false)
    success = putFromXML(coeff_ptr);
  else
      if(H5file!=true){
         APP_ABORT("Error in Opening HDF5");
      }
      else{
          hdf_archive hin(myComm);
    
          if(myComm->rank()==0){
            if (!hin.open(MOhref2)){
                APP_ABORT("SPOSetBase::putFromH5 missing or incorrect path to H5 file.");
            }
            //TO REVIEWERS:: IDEAL BEHAVIOUR SHOULD BE:
            /*
             if(!hin.push("PBC")
                 PBC=false;
             else
                if (!hin.read(PBC,"PBC"))
                    APP_ABORT("Could not read PBC dataset in H5 file. Probably corrupt file!!!.");
            // However, it always succeeds to enter the if condition even if the group does not exists...
            */             
            hin.push("PBC");
            PBC=false;
            hin.read(PBC,"PBC");
            hin.close();

          } 
          myComm->bcast(PBC);
          if (PBC)
             success = putPBCFromH5(MOhref2, coeff_ptr);
          else
             success = putMolFromH5(MOhref2, coeff_ptr);
      }
  bool success2 = transformSPOSet();
  if(debugc=="yes")
  {
    app_log() << "   Single-particle orbital coefficients dims="
      << C->rows() << " x " << C->cols() << std::endl;
    app_log() << C << std::endl;
  }

  init_LCOrbitalSetOpt(orbital_mix_magnitude);

  return success && success2;

 return true;
}

void SPOSetBase::checkObject()
{
  if(!(OrbitalSetSize == C->rows() && BasisSetSize == C->cols()))
  {
    app_error() << "   SPOSetBase::checkObject Linear coeffient for SPOSet is not consistent with the input." << std::endl;
    OHMMS::Controller->abort();
  }
}



bool SPOSetBase::putFromXML(xmlNodePtr coeff_ptr)
{
  Identity=true;
  int norbs=0;
  OhmmsAttributeSet aAttrib;
  aAttrib.add(norbs,"size");
  aAttrib.add(norbs,"orbitals");
  aAttrib.put(coeff_ptr);
  if(norbs < OrbitalSetSize)
  {
    return false;
    APP_ABORT("SPOSetBase::putFromXML missing or incorrect size");
  }
  if(norbs)
  {
    Identity=false;
    std::vector<ValueType> Ctemp;
    Ctemp.resize(norbs*BasisSetSize);
    setIdentity(Identity);
    putContent(Ctemp,coeff_ptr);
    int n=0,i=0;
    std::vector<ValueType>::iterator cit(Ctemp.begin());
    while(i<OrbitalSetSize)
    {
      if(Occ[n]>std::numeric_limits<RealType>::epsilon())
      {
        std::copy(cit,cit+BasisSetSize,(*C)[i]);
        i++;
      }
      n++;
      cit+=BasisSetSize;
    }
  }
  return true;
}

/** read data from a hdf5 file
 * @param norb number of orbitals to be initialized
 * @param fname hdf5 file name
 * @param occ_ptr xmlnode for occupation
 * @param coeff_ptr xmlnode for coefficients
 */
bool SPOSetBase::putMolFromH5(const char* fname, xmlNodePtr coeff_ptr)
{
#if defined(HAVE_LIBHDF5)
  ReportEngine PRE("SOA_SPOSetBase","SPOSetBase::putFromH5");
  int norbs=OrbitalSetSize;
  int neigs=BasisSetSize;
  int setVal=-1;
  std::string setname;
  OhmmsAttributeSet aAttrib;
  aAttrib.add(setVal,"spindataset");
  aAttrib.add(neigs,"size");
  aAttrib.add(neigs,"orbitals");
  aAttrib.put(coeff_ptr);
  setIdentity(false);
  hdf_archive hin(myComm);
  if(myComm->rank()==0){
    hin.open(fname);
    if (!hin.open(fname)){
        APP_ABORT("SPOSetBase::putFromH5 missing or incorrect path to H5 file.");
    }

    Matrix<RealType> Ctemp(BasisSetSize,BasisSetSize);
    char name[72];
    sprintf(name,"%s%d","/KPTS_0/eigenset_",setVal);
    setname=name;
    if(!hin.read(Ctemp,setname))
    {
       setname="SPOSetBase::putFromH5 Missing "+setname+" from HDF5 File.";
       APP_ABORT(setname.c_str());
    }
    hin.close();

    int n=0,i=0;
    while(i<norbs)
    {
      if(Occ[n]>0.0)
      {
        std::copy(Ctemp[n],Ctemp[n+1],(*C)[i]);
        i++;
      }
      n++;
    }
 }
 myComm->bcast(C->data(),C->size());
#else
  APP_ABORT("SPOSetBase::putFromH5 HDF5 is disabled.")
#endif
  return true;
}


/** read data from a hdf5 file
 * @param norb number of orbitals to be initialized
 * @param fname hdf5 file name
 * @param occ_ptr xmlnode for occupation
 * @param coeff_ptr xmlnode for coefficients
 */
bool SPOSetBase::putPBCFromH5(const char* fname, xmlNodePtr coeff_ptr)
{
#if defined(HAVE_LIBHDF5)
  ReportEngine PRE("SOA_SPOSetBase","SPOSetBase::putFromH5");
  int norbs=OrbitalSetSize;
  int neigs=BasisSetSize;
  int setVal=-1;
  int NbKpts;
  int KptIdx=0;
  bool IsComplex=false;
  PosType twist(0.0);
  PosType twistH5(0.0);
  std::string setname;
  OhmmsAttributeSet aAttrib;
  aAttrib.add(setVal,"spindataset");
  aAttrib.add(neigs,"size");
  aAttrib.add(neigs,"orbitals");
  aAttrib.put(coeff_ptr);
  setIdentity(false);
  hdf_archive hin(myComm);


  xmlNodePtr curtemp=coeff_ptr->parent->parent->parent;
  aAttrib.add(twist,"twist");
  aAttrib.put(curtemp);
  
  if(myComm->rank()==0){
    hin.open(fname);
    if (!hin.open(fname)){
        APP_ABORT("SPOSetBase::putFromH5 missing or incorrect path to H5 file.");
    }
    hin.push("parameters");
    hin.read(IsComplex,"IsComplex"); 
    hin.pop();
    hin.push("Nb_KPTS");
    hin.read(NbKpts,"Nbkpts"); 
    hin.pop();
    for (int i=0;i<NbKpts;i++)
    {
       char name[72];
       sprintf(name,"%s%d%s","/KPTS_",i,"/Coord");
       setname=name;
       hin.read(twistH5,setname);
       if(std::abs(twistH5[0]-twist[0])<1e-6 &&  std::abs(twistH5[1]-twist[1])<1e-6 && std::abs(twistH5[2]-twist[2])<1e-6)
       {
          KptIdx=i;
          break;
       }
    } 

    Matrix<RealType> Ctemp(BasisSetSize,BasisSetSize);

    char name[72];
    if(IsComplex){ 
         APP_ABORT("Complex Wavefunction not implemented yet. Please contact Developers");

    }
    else{
        sprintf(name,"%s%d%s%d","/KPTS_",KptIdx,"/eigenset_",setVal);
        setname=name;
        if(!hin.read(Ctemp,setname))
        {
           setname="SPOSetBase::putFromH5 Missing "+setname+" from HDF5 File.";
           APP_ABORT(setname.c_str());
        }
    }


#if defined (QMC_COMPLEX)
    APP_ABORT("Complex Wavefunction not implemented yet. Please contact Developers");
#endif //COMPLEX
    hin.close();

    int n=0,i=0;
    while(i<norbs)
    {
      if(Occ[n]>0.0)
      {
        std::copy(Ctemp[n],Ctemp[n+1],(*C)[i]);
        i++;
      }
      n++;
    }
 }
 myComm->bcast(C->data(),C->size());
#else
  APP_ABORT("SPOSetBase::putFromH5 HDF5 is disabled.")
#endif
  return true;
}


}
