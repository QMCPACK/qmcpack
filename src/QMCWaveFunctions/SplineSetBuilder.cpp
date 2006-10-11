//////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Jeongnim Kim
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
#include "OhmmsData/AttributeSet.h"
#include "Numerics/HDFSTLAttrib.h"
#include "QMCWaveFunctions/SplineSetBuilder.h"
//#include "Numerics/HDFTriCubicSpline.h"

namespace qmcplusplus {

  SplineSetBuilder::SplineSetBuilder(ParticleSet& p, PtclPoolType& psets):
    targetPtcl(p),ptclPool(psets) {
  }   

  /** Initialize the cubic grid
   */
  bool SplineSetBuilder::put(xmlNodePtr cur){
    cur=cur->children;
    while(cur != NULL) {
      string cname((const char*)(cur->name));
      if(cname == "spline") { //check this out
        if(GridXYZ==0) {
          GridXYZ = new GridType;
          GridXYZ->put(cur);
        }
      }
      cur=cur->next;
    }
    //grid
    return true;
  }

  /** create a SlaterDeterminant
   * @param cur xmlnode containing \<slaterdeterminant\>
   * @return a SlaterDeterminant
   *
   * @warning MultiSlaterDeterminant is not working yet.
   */
  SPOSetBase*
  SplineSetBuilder::createSPOSet(xmlNodePtr cur){

    string hrefname("NONE");
    int norb(0);
    OhmmsAttributeSet aAttrib;
    aAttrib.add(norb,"orbitals");
    aAttrib.add(hrefname,"href");
    aAttrib.put(cur);

    if(norb ==0) {
      app_error() << "SplineSetBuilder::createSPOSet failed. Check the attribte orbitals." << endl;
      return 0;
    }

    std::vector<int> npts(3);
    npts[0]=GridXYZ->nX;
    npts[1]=GridXYZ->nY;
    npts[2]=GridXYZ->nZ;
    std::vector<RealType> inData(npts[0]*npts[1]*npts[2]);

    SPOSetType* psi= new SPOSetType(norb);
    vector<int> occSet(norb);
    for(int i=0; i<norb; i++) occSet[i]=i;

    cur=cur->children;
    while(cur != NULL) {
      string cname((const char*)(cur->name));
      if(cname == "occupation") {
        string occ_mode("ground");
        const xmlChar* o=xmlGetProp(cur,(const xmlChar*)"mode");  
        if(o!= NULL) occ_mode = (const char*)o;
        //Do nothing if mode == ground
        if(occ_mode == "excited") {
          vector<int> occ_in, occRemoved;
          putContent(occ_in,cur);
          for(int k=0; k<occ_in.size(); k++) {
            if(occ_in[k]<0) 
              occRemoved.push_back(-occ_in[k]-1);
          }
          int kpopd=0;
          for(int k=0; k<occ_in.size(); k++) {
            if(occ_in[k]>0) 
              occSet[occRemoved[kpopd++]]=occ_in[k]-1;
          }
        }
        cout << "Index of the occupied states " << endl;
        for(int i=0; i< norb; i++) cout << occSet[i] << endl;
 
        hid_t h_file = H5Fopen(hrefname.c_str(),H5F_ACC_RDWR,H5P_DEFAULT);

        const xmlChar* h5path = xmlGetProp(cur,(const xmlChar*)"h5path");
        string hroot("/eigenstates_3/twist_0");
        if(h5path != NULL) hroot=(const char*)h5path;
        char wfname[128];
        for(int iorb=0; iorb<norb; iorb++) {
          sprintf(wfname,"%s/band_%d/eigenvector",hroot.c_str(),occSet[iorb]);
          app_log() << "   Reading spline function " << wfname << endl;
          SPOType* neworb=0;
          map<string,SPOType*>::iterator it(NumericalOrbitals.find(wfname));
          if(it == NumericalOrbitals.end()) {
            neworb=new SPOType(GridXYZ);
            HDFAttribIO<std::vector<RealType> > dummy(inData,npts);
            dummy.read(h_file,wfname);
            //////////////////////////////
            //HDFAttribIO<NumericalOrbitalType> dummy(*neworb);
            //dummy.read(h_file,"CubicSpline");
            ////////////////////////////
            neworb->reset(inData.begin(), inData.end(), targetPtcl.Lattice.BoxBConds[0]);
            NumericalOrbitals[wfname]=neworb;
          } else {
            neworb = (*it).second;
          }
          psi->add(neworb);
        }

        H5Fclose(h_file);
      }
      cur=cur->next;
    }
    return psi;
  }
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
