//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
/** @file PWOribitalBuilder.cpp
 * @brief Definition of a builder class for PWOrbitalSet
 */
#include "QMCWaveFunctions/PlaneWave/PWOrbitalBuilder.h"
#include "QMCWaveFunctions/PlaneWave/PWParameterSet.h"
#include "QMCWaveFunctions/Fermion/SlaterDet.h"
#include "QMCWaveFunctions/DummyBasisSet.h"
#include "QMCWaveFunctions/SPOSetScanner.h"
#include "OhmmsData/ParameterSet.h"
#include "OhmmsData/AttributeSet.h"
#include "Numerics/HDFSTLAttrib.h"
#include "Numerics/HDFNumericAttrib.h"
#include "Message/Communicate.h"

namespace qmcplusplus
{

PWOrbitalBuilder::PWOrbitalBuilder(ParticleSet& els, TrialWaveFunction& psi, PtclPoolType& psets)
  : OrbitalBuilderBase(els,psi), ptclPool(psets), hfileID(-1), rootNode(NULL)
#if !defined(EANBLE_SMARTPOINTER)
  ,myBasisSet(0)
#endif
{
  myParam=new PWParameterSet;
}

PWOrbitalBuilder::~PWOrbitalBuilder()
{
  delete myParam;
}

//All data parsing is handled here, outside storage classes.
bool PWOrbitalBuilder::put(xmlNodePtr cur)
{
  //save the parent
  rootNode=cur;
  //
  //Get wavefunction data and parameters from XML and HDF5
  //
  RealType ecut=-1.0;
  //close it if open
  if(hfileID>0)
    H5Fclose(hfileID);
  //check the current href
  hfileID = getH5(cur,"href");
  //no file, check the root
  if(hfileID<0)
    hfileID = getH5(rootNode,"href");
  bool success=true;
  //Move through the XML tree and read basis information
  cur = cur->children;
  while(cur != NULL)
  {
    std::string cname((const char*)(cur->name));
    if(cname == "basisset")
    {
      const xmlChar* aptr = xmlGetProp(cur,(const xmlChar*)"ecut");
      if(aptr != NULL)
        myParam->Ecut=atof((const char*)aptr);
    }
    else if(cname == "coefficients")
    {
      //close
      if(hfileID>0)
        H5Fclose(hfileID);
      hfileID=getH5(cur,"hdata");
    }
    else if(cname == sd_tag)
    {
      if(hfileID<0)
        hfileID=getH5(cur,"href");
      if(hfileID<0)
      {
        app_error() << "  Cannot create a SlaterDet due to missing h5 file" << std::endl;
        OHMMS::Controller->abort();
      }
      success = createPWBasis(cur);
      success = putSlaterDet(cur);
    }
    else if (cname == sposcanner_tag)
    {
      SPOSetScanner ascanner(spomap, targetPtcl, ptclPool);
      ascanner.put(cur);
    }
    cur=cur->next;
  }
  H5Fclose(hfileID);
  return success;
}


bool PWOrbitalBuilder::putSlaterDet(xmlNodePtr cur)
{
  //catch parameters
  myParam->put(cur);
  typedef SlaterDet SlaterDeterminant_t;
  typedef DiracDeterminantBase Det_t;
  SlaterDeterminant_t* sdet(new SlaterDeterminant_t(targetPtcl));
  int spin_group=0;
  cur=cur->children;
  while(cur != NULL)
  {
    std::string cname((const char*)(cur->name));
    //Which determinant?
    if(cname == "determinant")
    {
      std::string id("updet");
      std::string ref("0");
      OhmmsAttributeSet aAttrib;
      aAttrib.add(id,"id");
      aAttrib.add(ref,"ref");
      aAttrib.put(cur);
      if(ref == "0")
        ref=id;
      int firstIndex=targetPtcl.first(spin_group);
      std::map<std::string,SPOSetBasePtr>::iterator lit(spomap.find(ref));
      Det_t* adet=0;
      //int spin_group=0;
      if(lit == spomap.end())
      {
        app_log() << "  Create a PWOrbitalSet" << std::endl;;
        SPOSetBasePtr psi(createPW(cur,spin_group));
        sdet->add(psi,ref);
        spomap[ref] = psi;
        adet = new Det_t(psi,firstIndex);
      }
      else
      {
        app_log() << "  Reuse a PWOrbitalSet" << std::endl;
        adet = new Det_t((*lit).second,firstIndex);
      }
      app_log()<< "    spin=" << spin_group  << " id=" << id << " ref=" << ref << std::endl;
      if(adet)
      {
        adet->set(firstIndex,targetPtcl.last(spin_group)-firstIndex);
        sdet->add(adet,spin_group);
      }
      spin_group++;
    }
    cur = cur->next;
  }
  if(spin_group)
  {
    targetPsi.addOrbital(sdet,"SlaterDet");
  }
  else
  {
    APP_ABORT(" Failed to create a SlaterDet at PWOrbitalBuilder::putSlaterDet ");
  }
  return true;
}

/** The read routine - get data from XML and H5. Process it and build orbitals.
 *
 * - parameters
 *   -- num_tiwsts
 *   -- num_bands
 *   -- complex_coefficients
 *   -- maximum_ecut
 * - basis
 */
bool PWOrbitalBuilder::createPWBasis(xmlNodePtr cur)
{
  //recycle int and double reader
  int idata;
  double ddata;
  HDFAttribIO<int> hdfint(idata);
  HDFAttribIO<double> hdfdbl(ddata);
  //start of parameters
  hdfint.read(hfileID,"electrons/number_of_kpoints");
  int nkpts = idata;
  hdfint.read(hfileID,"electrons/number_of_spins");
  int nspin = idata;
  hdfint.read(hfileID,"electrons/kpoint_0/spin_0/number_of_states");
  int nbands = idata;
  myParam->numBands = nbands;
  app_log() << "Number of bands = " << nbands << std::endl;
  bool h5coefsreal = true;
  // Cutoff no longer present in the HDF file
  RealType h5ecut = 0.0;
  RealType ecut = 0.0;
  //end of parameters
  //check if input parameters are valid
  int nup=targetPtcl.last(0);
  int ndown=targetPtcl.getTotalNum()-nup;
  if(nbands < nup || nbands < ndown)
  {
    app_error() << "Not enough bands in h5 file" << std::endl;
    OHMMS::Controller->abort();
  }
  std::string tname=myParam->getTwistAngleName();
  TinyVector<double, OHMMS_DIM> TwistAngle_DP;
  HDFAttribIO<TinyVector<double, OHMMS_DIM> > hdfobj_twist(TwistAngle_DP);
  hdfobj_twist.read(hfileID,"/electrons/kpoint_0/reduced_k");
  TwistAngle=TwistAngle_DP;
#if defined(ENABLE_SMARTPOINTER)
  if(myBasisSet.get() ==0)
  {
    myBasisSet.reset(new PWBasis(TwistAngle));
  }
#else
  if(myBasisSet==0)
  {
    myBasisSet = new PWBasis(TwistAngle);
  }
#endif
  //Read the planewave basisset.
  //Note that the same data is opened here for each twist angle-avoids duplication in the
  //h5 file (which may become very large).
  //return the ecut to be used by the basis set
  RealType real_ecut = myParam->getEcut(ecut);
  //create at least one basis set but do resize the containers
  int nh5gvecs=myBasisSet->readbasis(hfileID,real_ecut,targetPtcl.Lattice,
                                     myParam->pwTag, myParam->pwMultTag);
  app_log() << "  num_twist = " << nkpts << std::endl;
  app_log() << "  twist angle = " << TwistAngle << std::endl;
  app_log() << "  num_bands = " << nbands << std::endl;
  app_log() << "  input maximum_ecut = " << ecut << std::endl;
  app_log() << "  current maximum_ecut = " << real_ecut << std::endl;
  app_log() << "  num_planewaves = " << nh5gvecs<< std::endl;
  return true;
}

SPOSetBase*
PWOrbitalBuilder::createPW(xmlNodePtr cur, int spinIndex)
{
  int nb=targetPtcl.last(spinIndex)-targetPtcl.first(spinIndex);
  std::vector<int> occBand(nb);
  for(int i=0; i<nb; i++)
    occBand[i]=i;
  typedef PWBasis::GIndex_t GIndex_t;
  GIndex_t nG(1);
  bool transform2grid=false;
  cur=cur->children;
  while(cur != NULL)
  {
    std::string cname((const char*)(cur->name));
    if(cname == "transform")
    {
      putContent(nG,cur);
      transform2grid=true;
    }
    else
      if(cname == "occupation")
      {
        std::string occMode("ground");
        int bandoffset(1);
        OhmmsAttributeSet aAttrib;
        aAttrib.add(spinIndex,"spindataset");
        aAttrib.add(occMode,"mode");
        aAttrib.add(bandoffset,"offset"); /* reserved for index offset */
        aAttrib.put(cur);
        if(occMode == "excited")
        {
          std::vector<int> occ;
          std::vector<int> deleted, added;
          putContent(occ,cur);
          for(int i=0; i<occ.size(); i++)
          {
            if(occ[i]<0)
              deleted.push_back(-occ[i]);
            else
              added.push_back(occ[i]);
          }
          if(deleted.size() != added.size())
          {
            app_error() << "  Numbers of deleted and added bands are not identical." << std::endl;
            OHMMS::Controller->abort();
          }
          for(int i=0; i<deleted.size(); i++)
          {
            occBand[deleted[i]-bandoffset]=added[i]-bandoffset;
          }
          app_log() << "  mode=\"excited\" Occupied states: " << std::endl;
          copy(occBand.begin(),occBand.end(),std::ostream_iterator<int>(app_log()," "));
          app_log() << std::endl;
        }
      }
    cur=cur->next;
  }
  //std::string tname=myParam->getTwistName();
  std::string tname="kpoint_0";
  //hid_t es_grp_id = H5Gopen(hfileID,myParam->eigTag.c_str());
  //hid_t es_grp_id = H5Gopen(hfileID,"electrons/kpoint_0/spin_0/state_0");
  hid_t es_grp_id = H5Gopen(hfileID,"electrons");
  
  //hid_t twist_grp_id = H5Gopen(es_grp_id,tname.c_str());
  hid_t twist_grp_id = H5Gopen(es_grp_id,"kpoint_0");
  //create a single-particle orbital set
  SPOSetType* psi=new SPOSetType;
  if(transform2grid)
  {
    nb=myParam->numBands;
    occBand.resize(nb);
    for(int i=0; i<nb; i++)
      occBand[i]=i;
  }
  //going to take care of occ
  psi->resize(myBasisSet,nb,true);
  if(myParam->hasComplexData(hfileID))//input is complex
  {
    //app_log() << "  PW coefficients are complex." << std::endl;
    typedef std::vector<std::complex<RealType> > TempVecType;
    typedef std::vector<std::complex<double> > TempVecType_DP;
    TempVecType_DP coefs_DP(myBasisSet->inputmap.size());
    HDFAttribIO<TempVecType_DP> hdfobj_coefs(coefs_DP);
    int ib=0;
    while(ib<nb)
    {
      std::string bname(myParam->getBandName(occBand[ib],spinIndex));
      app_log() << "  Reading " << tname <<"/"<< bname << std::endl;
      hid_t band_grp_id =  H5Gopen(twist_grp_id,bname.c_str());
      hdfobj_coefs.read(band_grp_id,"psi_g");
      TempVecType coefs(coefs_DP.begin(), coefs_DP.end());
      psi->addVector(coefs,ib);
      H5Gclose(band_grp_id);
      ++ib;
    }
  }
  else
  {
    // It appears the coefficients are always stored as complex in the HDF file?
    //app_log() << "  PW coefficients are real." << std::endl;
    typedef std::vector<std::complex<RealType> > ComplexTempVecType;
    typedef std::vector<std::complex<double> > ComplexTempVecType_DP;
    ComplexTempVecType_DP complex_coefs_DP(myBasisSet->inputmap.size());
    HDFAttribIO<ComplexTempVecType_DP> hdfobj_complex_coefs(complex_coefs_DP);
    int ib=0;
    while(ib<nb)
    {
      std::string bname(myParam->getBandName(occBand[ib],spinIndex));
      app_log() << "  Reading " << tname <<"/"<< bname << std::endl;
      hid_t band_grp_id =  H5Gopen(twist_grp_id,bname.c_str());
      hdfobj_complex_coefs.read(band_grp_id,"psi_g");
      ComplexTempVecType complex_coefs(complex_coefs_DP.begin(),complex_coefs_DP.end());
      psi->addVector(complex_coefs,ib);
      H5Gclose(band_grp_id);
      ++ib;
    }
  }
  H5Gclose(twist_grp_id);
  H5Gclose(es_grp_id);
#if defined(QMC_COMPLEX)
  if(transform2grid)
  {
    app_warning() << "  Going to transform on grid " << std::endl;
    transform2GridData(nG,spinIndex,*psi);
  }
#endif
  return psi;
}

#if defined(QMC_COMPLEX)
void PWOrbitalBuilder::transform2GridData(PWBasis::GIndex_t& nG, int spinIndex, PWOrbitalSet& pwFunc)
{
  std::ostringstream splineTag;
  splineTag << "eigenstates_"<<nG[0]<<"_"<<nG[1]<<"_"<<nG[2];
  herr_t status = H5Eset_auto(NULL, NULL);
  app_log() << " splineTag " << splineTag.str() << std::endl;
  hid_t es_grp_id;
  status = H5Gget_objinfo (hfileID, splineTag.str().c_str(), 0, NULL);
  if(status)
  {
    es_grp_id = H5Gcreate(hfileID,splineTag.str().c_str(),0);
    HDFAttribIO<PWBasis::GIndex_t> t(nG);
    t.write(es_grp_id,"grid");
  }
  else
  {
    es_grp_id = H5Gopen(hfileID,splineTag.str().c_str());
  }
  std::string tname=myParam->getTwistName();
  hid_t twist_grp_id;
  status = H5Gget_objinfo (es_grp_id, tname.c_str(), 0, NULL);
  if(status)
    twist_grp_id = H5Gcreate(es_grp_id,tname.c_str(),0);
  else
    twist_grp_id = H5Gopen(es_grp_id,tname.c_str());
  TinyVector<double, OHMMS_DIM> TwistAngle_DP;
  TwistAngle_DP=TwistAngle;
  HDFAttribIO<TinyVector<double, OHMMS_DIM> > hdfobj_twist(TwistAngle_DP);
  hdfobj_twist.write(twist_grp_id,"twist_angle");
  ParticleSet::ParticleLayout_t& lattice(targetPtcl.Lattice);
  RealType dx=1.0/static_cast<RealType>(nG[0]-1);
  RealType dy=1.0/static_cast<RealType>(nG[1]-1);
  RealType dz=1.0/static_cast<RealType>(nG[2]-1);
#if defined(VERYTINYMEMORY)
  typedef Array<ParticleSet::ParticleValue_t,3> StorageType;
  StorageType inData(nG[0],nG[1],nG[2]);
  int ib=0;
  while(ib<myParam->numBands)
  {
    std::string bname(myParam->getBandName(ib));
    status = H5Gget_objinfo (twist_grp_id, bname.c_str(), 0, NULL);
    hid_t band_grp_id, spin_grp_id=-1;
    if(status)
    {
      band_grp_id =  H5Gcreate(twist_grp_id,bname.c_str(),0);
    }
    else
    {
      band_grp_id =  H5Gopen(twist_grp_id,bname.c_str());
    }
    hid_t parent_id=band_grp_id;
    if(myParam->hasSpin)
    {
      bname=myParam->getSpinName(spinIndex);
      status = H5Gget_objinfo (band_grp_id, bname.c_str(), 0, NULL);
      if(status)
      {
        spin_grp_id =  H5Gcreate(band_grp_id,bname.c_str(),0);
      }
      else
      {
        spin_grp_id =  H5Gopen(band_grp_id,bname.c_str());
      }
      parent_id=spin_grp_id;
    }
    for(int ig=0; ig<nG[0]; ig++)
    {
      RealType x=ig*dx;
      for(int jg=0; jg<nG[1]; jg++)
      {
        RealType y=jg*dy;
        for(int kg=0; kg<nG[2]; kg++)
        {
          inData(ig,jg,kg)= pwFunc.evaluate(ib,lattice.toCart(PosType(x,y,kg*dz)));
        }
      }
    }
    app_log() << "  Add spline data " << ib << " h5path=" << tname << "/eigvector" << std::endl;
    HDFAttribIO<StorageType> t(inData);
    t.write(parent_id,myParam->eigvecTag.c_str());
    if(spin_grp_id>=0)
      H5Gclose(spin_grp_id);
    H5Gclose(band_grp_id);
    ++ib;
  }
#else
  typedef Array<ParticleSet::ParticleValue_t,3> StorageType;
  std::vector<StorageType*> inData;
  int nb=myParam->numBands;
  for(int ib=0; ib<nb; ib++)
    inData.push_back(new StorageType(nG[0],nG[1],nG[2]));
  PosType tAngle=targetPtcl.Lattice.k_cart(TwistAngle);
  PWOrbitalSet::ValueVector_t phi(nb);
  for(int ig=0; ig<nG[0]; ig++)
  {
    RealType x=ig*dx;
    for(int jg=0; jg<nG[1]; jg++)
    {
      RealType y=jg*dy;
      for(int kg=0; kg<nG[2]; kg++)
      {
        targetPtcl.R[0]=lattice.toCart(PosType(x,y,kg*dz));
        pwFunc.evaluate(targetPtcl,0,phi);
        RealType x(dot(targetPtcl.R[0],tAngle));
        ValueType phase(std::cos(x),-std::sin(x));
        for(int ib=0; ib<nb; ib++)
          (*inData[ib])(ig,jg,kg)=phase*phi[ib];
      }
    }
  }
  for(int ib=0; ib<nb; ib++)
  {
    std::string bname(myParam->getBandName(ib));
    status = H5Gget_objinfo (twist_grp_id, bname.c_str(), 0, NULL);
    hid_t band_grp_id, spin_grp_id=-1;
    if(status)
    {
      band_grp_id =  H5Gcreate(twist_grp_id,bname.c_str(),0);
    }
    else
    {
      band_grp_id =  H5Gopen(twist_grp_id,bname.c_str());
    }
    hid_t parent_id=band_grp_id;
    if(myParam->hasSpin)
    {
      bname=myParam->getSpinName(spinIndex);
      status = H5Gget_objinfo (band_grp_id, bname.c_str(), 0, NULL);
      if(status)
      {
        spin_grp_id =  H5Gcreate(band_grp_id,bname.c_str(),0);
      }
      else
      {
        spin_grp_id =  H5Gopen(band_grp_id,bname.c_str());
      }
      parent_id=spin_grp_id;
    }
    app_log() << "  Add spline data " << ib << " h5path=" << tname << "/eigvector" << std::endl;
    HDFAttribIO<StorageType> t(*(inData[ib]));
    t.write(parent_id,myParam->eigvecTag.c_str());
    if(spin_grp_id>=0)
      H5Gclose(spin_grp_id);
    H5Gclose(band_grp_id);
  }
  for(int ib=0; ib<nb; ib++)
    delete inData[ib];
#endif
  H5Gclose(twist_grp_id);
  H5Gclose(es_grp_id);
}
#endif

hid_t PWOrbitalBuilder::getH5(xmlNodePtr cur, const char* aname)
{
  const xmlChar* aptr = xmlGetProp(cur,(const xmlChar*)aname);
  if(aptr == NULL)
  {
    return -1;
  }
  hid_t h = H5Fopen((const char*)aptr,H5F_ACC_RDWR,H5P_DEFAULT);
  if(h<0)
  {
    app_error() << " Cannot open " << (const char*)aptr << " file." << std::endl;
    OHMMS::Controller->abort();
  }
  myParam->checkVersion(h);
  //overwrite the parameters
  myParam->put(rootNode);
  return h;
}

}
