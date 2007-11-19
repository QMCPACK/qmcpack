//////////////////////////////////////////////////////////////////
// (c) Copyright 2006-  by Jeongnim Kim
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
#include "QMCWaveFunctions/SPOSetBase.h"
#if defined(HAVE_LIBHDF5)
#include "Numerics/HDFNumericAttrib.h"
#endif
#include "OhmmsData/AttributeSet.h"
#include "Message/Communicate.h"
#include <limits>

namespace qmcplusplus {

  /** Parse the xml file for information on the Dirac determinants.
   *@param cur the current xmlNode
   */
  bool SPOSetBase::put(xmlNodePtr cur) {
    //initialize the number of orbital by the basis set size
    int norb= BasisSetSize;

    OhmmsAttributeSet aAttrib;
    aAttrib.add(norb,"orbitals"); aAttrib.add(norb,"size");
    aAttrib.put(cur);
    //const xmlChar* norb_ptr=xmlGetProp(cur, (const xmlChar *)"orbitals");
    //if(norb_ptr != NULL) { 
    //  norb=atoi((const char*)norb_ptr);
    //}

    setOrbitalSetSize(norb);

    const xmlChar* h=xmlGetProp(cur, (const xmlChar*)"href");
    xmlNodePtr occ_ptr=NULL;
    xmlNodePtr coeff_ptr=NULL;
    cur = cur->xmlChildrenNode;
    while(cur != NULL) {
      string cname((const char*)(cur->name));
      if(cname == "occupation") {
        occ_ptr=cur;
      } else if(cname.find("coeff") < cname.size() || cname == "parameter" || cname == "Var") {
        coeff_ptr=cur;
      }
      cur=cur->next;
    }

    if(coeff_ptr == NULL) {
      app_log() << "   Using Identity for the LCOrbitalSet " << endl;
      return setIdentity(true);
    }

    bool success=putOccupation(occ_ptr);

    if(h == NULL) {
      success = putFromXML(coeff_ptr);
    } else {
      success = putFromH5((const char*)h, coeff_ptr);
    }
    return success;
  }

  void SPOSetBase::checkObject() {
    if(!(OrbitalSetSize == C.rows() && BasisSetSize == C.cols()))
    {
      app_error() << "   SPOSetBase::checkObject Linear coeffient for SPOSet is not consistent with the input." << endl; 
      OHMMS::Controller->abort();
    }
  }

  bool SPOSetBase::putOccupation(xmlNodePtr occ_ptr) {

    //die??
    if(BasisSetSize ==0) return false;

    Occ.resize(BasisSetSize,0.0);
    for(int i=0; i<OrbitalSetSize; i++) Occ[i]=1.0;

    vector<int> occ_in;
    string occ_mode("table");
    if(occ_ptr == NULL) {
      occ_mode="ground";
    } else {
      const xmlChar* o=xmlGetProp(occ_ptr,(const xmlChar*)"mode");  
      if(o) occ_mode = (const char*)o;
    }
    //Do nothing if mode == ground
    if(occ_mode == "excited") {
      putContent(occ_in,occ_ptr);
      for(int k=0; k<occ_in.size(); k++) {
        if(occ_in[k]<0) //remove this, -1 is to adjust the base
          Occ[-occ_in[k]-1]=0.0;
        else
          Occ[occ_in[k]-1]=1.0;
      }
    } else if(occ_mode == "table") {
      putContent(Occ,occ_ptr);
    }

    return true;
  }

  bool SPOSetBase::putFromXML(xmlNodePtr coeff_ptr) {
    vector<ValueType> Ctemp;
    int total(OrbitalSetSize);
    Identity=true;
    if(coeff_ptr != NULL){
      if(xmlHasProp(coeff_ptr, (const xmlChar*)"size")){
        total = atoi((const char*)(xmlGetProp(coeff_ptr, (const xmlChar *)"size")));
      } 
      Ctemp.resize(total*BasisSetSize);
      putContent(Ctemp,coeff_ptr);
      Identity = false;
    }

    setIdentity(Identity);

    if(!Identity) {
      int n=0,i=0;
      vector<ValueType>::iterator cit(Ctemp.begin());
      while(i<OrbitalSetSize){
        if(Occ[n]>numeric_limits<RealType>::epsilon()){
          std::copy(cit,cit+BasisSetSize,C[i]);
          i++; 
        }
        n++;cit+=BasisSetSize;
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
    bool SPOSetBase::putFromH5(const char* fname, xmlNodePtr coeff_ptr) {

#if defined(HAVE_LIBHDF5)
      int norbs=OrbitalSetSize;
      int neigs=BasisSetSize;

      string setname("invalid");
      if(coeff_ptr) {
         const xmlChar* d=xmlGetProp(coeff_ptr,(const xmlChar*)"dataset");  
         if(d) setname = (const char*)d;
         d=xmlGetProp(coeff_ptr,(const xmlChar*)"size");  
         if(d) neigs=atoi((const char*)d);
      }

      setIdentity(false);

      if(setname != "invalid") {
        Matrix<RealType> Ctemp(neigs,BasisSetSize);
        hid_t h_file=  H5Fopen(fname,H5F_ACC_RDWR,H5P_DEFAULT);
        HDFAttribIO<Matrix<RealType> > h(Ctemp);
        h.read(h_file,setname.c_str());
        H5Fclose(h_file);
	int n=0,i=0;
	while(i<norbs){
	  if(Occ[n]>numeric_limits<RealType>::epsilon()){
            std::copy(Ctemp[n],Ctemp[n+1],C[i]);
	    i++; 
	  }
	  n++;
	}
      }

#else
      ERRORMSG("HDF5 is disabled. Using identity")
#endif
      return true;
    }

}

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

