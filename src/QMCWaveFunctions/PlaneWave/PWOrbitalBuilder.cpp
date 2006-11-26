/////////////////////////////////////////////////////////////////
// (c) Copyright 2006-  by Kris Delaney and Jeongnim Kim
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
/** @file PWOribitalBuilder.cpp
 * @brief Definition of a builder class for PWOrbitalSet
 */
#include "Utilities/OhmmsInfo.h"
#include "QMCWaveFunctions/PlaneWave/PWOrbitalBuilder.h"
#include "QMCWaveFunctions/PlaneWave/PWParameterSet.h"
#include "QMCWaveFunctions/Fermion/SlaterDet.h"
#include "QMCWaveFunctions/DummyBasisSet.h"
#include "OhmmsData/ParameterSet.h"
#include "OhmmsData/AttributeSet.h"
#include "Numerics/HDFSTLAttrib.h"
#include "Numerics/HDFNumericAttrib.h"
#include "Message/Communicate.h"

namespace qmcplusplus {

  PWOrbitalBuilder::PWOrbitalBuilder(ParticleSet& els, TrialWaveFunction& psi):
    OrbitalBuilderBase(els,psi),hfileID(-1)
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
  bool PWOrbitalBuilder::put(xmlNodePtr cur) {

    //catch the parameters at the root level
    myParam->put(cur);
    //
    //Get wavefunction data and parameters from XML and HDF5
    //
    RealType ecut=-1.0;
    
    if(hfileID>0) H5Fclose(hfileID);

    hfileID = getH5(cur,"href");

    bool success=true;
    //Move through the XML tree and read basis information
    cur = cur->children;
    while(cur != NULL) {
      string cname((const char*)(cur->name));
      if(cname == "basisset") 
      {
        const xmlChar* aptr = xmlGetProp(cur,(const xmlChar*)"ecut");
        if(aptr != NULL) ecut=atof((const char*)aptr);
        myParam->put(cur);
      } 
      else if(cname == "coefficients")
      {
        hfileID=getH5(cur,"hdata");
      } 
      else if(cname == "slaterdeterminant") {
        if(hfileID<0)
          hfileID=getH5(cur,"href");

        if(hfileID<0)
        {
          app_error() << "  Cannot create a SlaterDet due to missing h5 file" << endl;
          OHMMS::Controller->abort();
        }

        success = createPWBasis(cur);
        success = putSlaterDet(cur);
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

    SlaterDeterminant_t* sdet(new SlaterDeterminant_t);

    int sid=0;
    int firstIndex=0;
    cur=cur->children;
    while(cur != NULL) {
      string cname((const char*)(cur->name));
      //Which determinant?
      if(cname == "determinant") {
        string id("updet");
        string ref("0");
        OhmmsAttributeSet aAttrib;
        aAttrib.add(id,"id");
        aAttrib.add(ref,"ref");
        aAttrib.put(cur);

        if(ref == "0") ref=id;

        map<string,SPOSetBasePtr>::iterator lit(PWOSet.find(ref));
        Det_t* adet=0;
        int nstates=0;
        if(lit == PWOSet.end())
        {
          app_log() << "  Create a PWOrbitalSet" << endl;;
          SPOSetBasePtr psi(createPW(cur,sid));
          PWOSet[id]=psi;
          adet= new Det_t(psi,firstIndex);
          nstates=psi->getOrbitalSetSize();
        }
        else
        {
          app_log() << "  Reuse a PWOrbitalSet" << endl;
          adet= new Det_t((*lit).second,firstIndex);
          nstates=(*lit).second->getOrbitalSetSize();
        }

        app_log()<< "    spin=" << sid  << " id=" << id << " ref=" << ref << endl; 

        if(adet) 
        {
          adet->set(firstIndex,nstates);
          sdet->add(adet);
        }

        firstIndex += nstates;
        ++sid;
      }
      cur = cur->next;
    }


    if(sid)
    {
      targetPsi.addOrbital(sdet);
    }
    else
    {
      app_error() << "  Failed to create a SlaterDet at PWOrbitalBuilder::putSlaterDet " << endl;
      OHMMS::Controller->abort();
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
    hid_t grp_id = H5Gopen(hfileID,myParam->paramTag.c_str());
    
    hdfint.read(grp_id,"num_twists");
    int nkpts = idata; 

    hdfint.read(grp_id,"num_bands");
    int nbands = idata;
    myParam->numBands = nbands;

    hdfint.read(grp_id,"complex_coefficients");
    bool h5coefsreal = (idata==0);

    hdfdbl.read(grp_id,"maximum_ecut");
    RealType h5ecut = ddata;
    RealType ecut = ddata;
    
    H5Gclose(grp_id);
    //end of parameters

    //check if input parameters are valid
    int nup=targetPtcl.last(0);
    int ndown=targetPtcl.getTotalNum()-nup;
    if(nbands < nup || nbands < ndown){
      app_error() << "Not enough bands in h5 file" << endl;
      OHMMS::Controller->abort();
    }

    string tname=myParam->getTwistAngleName();
    //hid_t es_grp_id = H5Gopen(hfile,myParam->eigTag.c_str());
    //hid_t twist_grp_id = H5Gopen(es_grp_id,tname.c_str());
    PosType twist;
    HDFAttribIO<PosType> hdfobj_twist(twist);
    hdfobj_twist.read(hfileID,tname.c_str());

#if defined(ENABLE_SMARTPOINTER)
    if(myBasisSet.get() ==0)
    {
      myBasisSet.reset(new PWBasis(twist));
    }
#else
    if(myBasisSet==0) {
      myBasisSet = new PWBasis(twist);
    }
#endif
    //Read the planewave basisset.
    //Note that the same data is opened here for each twist angle-avoids duplication in the
    //h5 file (which may become very large).

    grp_id = H5Gopen(hfileID,myParam->basisTag.c_str());
    //create at least one basis set but do resize the containers
    int nh5gvecs=myBasisSet->readbasis(grp_id,ecut,targetPtcl.Lattice,myParam->pwTag);
    H5Gclose(grp_id); //Close PW Basis group

    app_log() << "  num_twist = " << nkpts << endl;
    app_log() << "  twist angle = " << twist << endl;
    app_log() << "  num_bands = " << nbands << endl;
    app_log() << "  maximum_ecut = " << ecut << endl;
    app_log() << "  num_planewaves = " << nh5gvecs<< endl;

    return true;
  }

  SPOSetBase* 
    PWOrbitalBuilder::createPW(xmlNodePtr cur, int spinIndex)
  {

    int nb=targetPtcl.last(spinIndex)-targetPtcl.first(spinIndex);
    vector<int> occBand(nb);
    for(int i=0;i<nb; i++) occBand[i]=i;

    typedef PWBasis::GIndex_t GIndex_t;
    GIndex_t nG(1);
    bool transform2grid=false;

    cur=cur->children;
    while(cur != NULL)
    {
      string cname((const char*)(cur->name));
      if(cname == "transform")
      {
        putContent(nG,cur);
        transform2grid=true;
      }
      else if(cname == "occupation")
      {
        string occMode("ground");
        int bandoffset(1);
        OhmmsAttributeSet aAttrib;
        aAttrib.add(spinIndex,"spindataset");
        aAttrib.add(occMode,"mode");
        aAttrib.add(bandoffset,"offset"); /* reserved for index offset */
        aAttrib.put(cur);

        if(occMode == "excited")
        {
          vector<int> occ;
          vector<int> deleted, added;

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
            app_error() << "  Numbers of deleted and added bands are not identical." << endl;
            OHMMS::Controller->abort();
          }

          for(int i=0; i<deleted.size(); i++)
          {
            occBand[deleted[i]-bandoffset]=added[i]-bandoffset;
          }

          app_log() << "  mode=\"excited\" Occupied states: " << endl;
          std::copy(occBand.begin(),occBand.end(),ostream_iterator<int>(app_log()," "));
          app_log() << endl;
        }
      } 
      cur=cur->next;
    }

    string tname=myParam->getTwistName();
    hid_t es_grp_id = H5Gopen(hfileID,myParam->eigTag.c_str());
    hid_t twist_grp_id = H5Gopen(es_grp_id,tname.c_str());

    ///create PWOrbital
    PWOrbitalSet* psi=new PWOrbitalSet;

    if(transform2grid)
    {
      nb=myParam->numBands;
      occBand.resize(nb);
      for(int i=0;i<nb; i++) occBand[i]=i;
    }

    //going to take care of occ
    psi->resize(myBasisSet,nb,true);

    typedef std::vector<ValueType> TempVecType;
    TempVecType coefs;
    HDFAttribIO<TempVecType> hdfobj_coefs(coefs);

    int ib=0;
    while(ib<nb) 
    {
      string bname(myParam->getBandName(occBand[ib],spinIndex));
      hid_t band_grp_id =  H5Gopen(twist_grp_id,bname.c_str());
      hdfobj_coefs.read(band_grp_id,myParam->eigvecTag.c_str());
      psi->addVector(coefs,ib);
      H5Gclose(band_grp_id);
      ++ib;
    }

    H5Gclose(twist_grp_id);
    H5Gclose(es_grp_id);

    if(transform2grid)
    {
      app_warning() << "  Going to transform on grid " << endl;
      transform2GridData(nG,spinIndex,*psi);
    }

    return psi;
  }

  void PWOrbitalBuilder::transform2GridData(PWBasis::GIndex_t& nG, int spinIndex, PWOrbitalSet& pwFunc)
  {
    ostringstream splineTag;
    splineTag << "eigenstates_"<<nG[0]<<"_"<<nG[1]<<"_"<<nG[2];
    herr_t status = H5Eset_auto(NULL, NULL);

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

    string tname=myParam->getTwistName();
    hid_t twist_grp_id;
    status = H5Gget_objinfo (es_grp_id, tname.c_str(), 0, NULL);
    if(status)
      twist_grp_id = H5Gcreate(es_grp_id,tname.c_str(),0);
    else
      twist_grp_id = H5Gopen(es_grp_id,tname.c_str());

    int ib=0;
    while(ib<myParam->numBands) 
    {
      string bname(myParam->getBandName(ib));
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

      HDFAttribIO<PWBasis::GIndex_t> t(nG);
      t.write(parent_id,myParam->eigvecTag.c_str());

      if(spin_grp_id>=0) H5Gclose(spin_grp_id);
      H5Gclose(band_grp_id);
      ++ib;
    }

    H5Gclose(twist_grp_id);
    H5Gclose(es_grp_id);
  }

  hid_t PWOrbitalBuilder::getH5(xmlNodePtr cur, const char* aname)
  {
    const xmlChar* aptr = xmlGetProp(cur,(const xmlChar*)aname);
    if(aptr == NULL){
      return -1;
    }

    hid_t h = H5Fopen((const char*)aptr,H5F_ACC_RDWR,H5P_DEFAULT);
    if(h<0)
    {
      app_error() << " Cannot open " << (const char*)aptr << " file." << endl;
      OHMMS::Controller->abort();
    }

    myParam->checkVersion(h);

    return h;
  }
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/
