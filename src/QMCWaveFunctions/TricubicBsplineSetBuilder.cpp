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
#include "QMCWaveFunctions/TricubicBsplineSetBuilder.h"

namespace qmcplusplus {

  TricubicBsplineSetBuilder::TricubicBsplineSetBuilder(ParticleSet& p, PtclPoolType& psets):
    targetPtcl(p),ptclPool(psets),GridXYZ(0) {
  }   

  /** Initialize the cubic grid
   */
  bool TricubicBsplineSetBuilder::put(xmlNodePtr cur){
    if(GridXYZ==0) {
      GridXYZ = new GridType;
      GridXYZ->put(cur);
    }
    return true;
  }

  /** Specialization for Array<double,3> */
  template<>
    struct HDFAttribIO<Array<double,3> >: public HDFAttribIOBase 
    {
      typedef Array<double,3> ArrayType_t;
      ArrayType_t&  ref;
      HDFAttribIO<ArrayType_t>(ArrayType_t& a):ref(a) { }
      inline void write(hid_t grp, const char* name) 
      {
      }

      inline void read(hid_t grp, const char* name) 
      {
        std::vector<int> npts(3);
        npts[0]=ref.size(0);
        npts[1]=ref.size(1);
        npts[2]=ref.size(2);
        HDFAttribIO<std::vector<double> > dummy(ref.storage(),npts);
        dummy.read(grp,name);
      }
    };

  /** create a SlaterDeterminant
   * @param cur xmlnode containing \<slaterdeterminant\>
   * @return a SlaterDeterminant
   *
   * @warning MultiSlaterDeterminant is not working yet.
   */
  SPOSetBase*
  TricubicBsplineSetBuilder::createSPOSet(xmlNodePtr cur){

    string hrefname("NONE");
    int norb(0);
    int degeneracy(1);
    OhmmsAttributeSet aAttrib;
    aAttrib.add(norb,"orbitals");
    aAttrib.add(degeneracy,"degeneracy");
    aAttrib.add(hrefname,"href");
    aAttrib.put(cur);

    if(norb ==0) {
      app_error() << "TricubicBsplineSetBuilder::createSPOSet failed. Check the attribte orbitals." << endl;
      return 0;
    }

    app_log() << "    Only valid for spin-unpolarized cases " << endl;
    app_log() << "    Degeneracy = " << degeneracy << endl;
    std::vector<int> npts(3);
    npts[0]=GridXYZ->nX;
    npts[1]=GridXYZ->nY;
    npts[2]=GridXYZ->nZ;
    StorageType inData(npts[0],npts[1],npts[2]);

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


        const xmlChar* h5path = xmlGetProp(cur,(const xmlChar*)"h5path");
        string hroot("/eigenstates_3/twist_0");
        if(h5path != NULL) hroot=(const char*)h5path;
        char wfname[128],wfshortname[16];
        
        map<string,OrbitalGroupType*>::iterator git(myBasis.find("0"));
        OrbitalGroupType *abasis=0;
        if(git == myBasis.end())
        {
          abasis = new OrbitalGroupType;
          abasis->setGrid(GridXYZ->x_min,GridXYZ->x_max,
              GridXYZ->y_min,GridXYZ->y_max, GridXYZ->z_min,GridXYZ->z_max,
              GridXYZ->nX-1, GridXYZ->nY-1, GridXYZ->nZ-1);
          myBasis["0"]=abasis;
        } 
        else
        {
          abasis=(*git).second;
        }

        hid_t h_file = H5Fopen(hrefname.c_str(),H5F_ACC_RDWR,H5P_DEFAULT);
        for(int iorb=0; iorb<norb; iorb++) {
          sprintf(wfname,"%s/band_%d/eigenvector",hroot.c_str(),occSet[iorb]/degeneracy);
          sprintf(wfshortname,"b%d",occSet[iorb]/degeneracy);
          StorageType* newP=0;
          map<string,StorageType*>::iterator it(BigDataSet.find(wfshortname));
          if(it == BigDataSet.end()) {
            newP=new StorageType;
            HDFAttribIO<StorageType> dummy(inData);
            dummy.read(h_file,wfname);
            BigDataSet[wfshortname]=newP;
            abasis->add(iorb,inData,newP);
            app_log() << "   Reading spline function " << wfname << endl;
            //PosType r0(4.9194332197e+00,4.5695928280e+00,1.2260692483e+01);
            //double val=neworb->evaluate(r0);
            //abort();
          } else {
            app_log() << "   Reusing spline function " << wfname << endl;
          }
          //psi->add(neworb);
        }
        psi->add(abasis);
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
