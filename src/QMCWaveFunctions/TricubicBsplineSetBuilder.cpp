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
#include "QMCWaveFunctions/PlaneWave/PWParameterSet.h"
#include "QMCWaveFunctions/TricubicBsplineSetBuilder.h"
#include "Numerics/TricubicBsplineSet.h"
#include "Numerics/TricubicBsplineGCSet.h"
#include "QMCWaveFunctions/GroupedOrbitalSet.h"
#include "QMCWaveFunctions/ElectronGas/ElectronGasOrbitalBuilder.h"
#include "QMCWaveFunctions/ElectronGas/HEGGrid.h"
#include "Numerics/OhmmsBlas.h"
#include "Message/Communicate.h"
//#define DEBUG_BSPLINE_EG
namespace qmcplusplus {

  //initialize the static data member
  map<string,TricubicBsplineSetBuilder::StorageType*> TricubicBsplineSetBuilder::BigDataSet;

  TricubicBsplineSetBuilder::TricubicBsplineSetBuilder(ParticleSet& p, PtclPoolType& psets, xmlNodePtr cur):
    targetPtcl(p),ptclPool(psets),rootNode(cur), LowerBox(0.0),UpperBox(1.0),BoxGrid(2)
    {
      for(int idim=0; idim<DIM; idim++)
        UpperBox[idim]=targetPtcl.Lattice.R(idim,idim);
      myParam=new PWParameterSet;
  }   

  TricubicBsplineSetBuilder::~TricubicBsplineSetBuilder()
  {
    delete myParam;
  }

  /** Initialize the cubic grid
   */
  bool TricubicBsplineSetBuilder::put(xmlNodePtr cur){

    curH5Fname.clear();

    //save the name of current HDF5 file name
    OhmmsAttributeSet aAttrib;
    aAttrib.add(curH5Fname,"href");
    aAttrib.put(rootNode);

    OpenEndGrid=false;
    cur=cur->children;
    while(cur != NULL)
    {
      std::string cname((const char*)(cur->name));
      if(cname == "grid") {
        string closedEnd("yes");
        int idir=0,npts=2; 
        RealType ri=0,rf=-1.0;
        OhmmsAttributeSet aAttrib;
        aAttrib.add(idir,"dir");
        aAttrib.add(ri,"ri");
        aAttrib.add(rf,"rf");
        aAttrib.add(npts,"npts");
        aAttrib.add(closedEnd,"closed");
        aAttrib.put(cur);
        //bound check
        RealType del=rf-ri;
        if(del>0 && del<UpperBox[idir]) UpperBox[idir]=rf;
        LowerBox[idir]=ri;
        BoxGrid[idir]=npts;
        OpenEndGrid = (closedEnd != "yes");
      }
      cur=cur->next;
    }
    app_log() << "  TricubicBsplineSetBuilder set a global box " << endl;
    app_log() << "    Lower box " << LowerBox << endl;
    app_log() << "    Upper box " << UpperBox << endl;
    app_log() << "    Box grid "  << BoxGrid << endl;
    return true;
  }

  template<typename OGT>
    SPOSetBase* TricubicBsplineSetBuilder::createBsplineBasisSet(xmlNodePtr cur, OGT* abasis)
  {
    int norb(0);
    int degeneracy(1);
    OhmmsAttributeSet aAttrib;
    aAttrib.add(norb,"orbitals"); aAttrib.add(norb,"size");
    aAttrib.add(degeneracy,"degeneracy");
    aAttrib.put(cur);

    vector<int> occSet(norb);
    for(int i=0; i<norb; i++) occSet[i]=i;

    //set the root name
    char hroot[128];
    sprintf(hroot,"/%s/%s%d",myParam->eigTag.c_str(),myParam->twistTag.c_str(),myParam->twistIndex);

    int spinIndex=0;
    cur=cur->children;
    while(cur != NULL) {
      string cname((const char*)(cur->name));
      if(cname == "occupation") {
        string occ_mode("ground");
        OhmmsAttributeSet oAttrib;
        oAttrib.add(occ_mode,"mode");
        oAttrib.add(spinIndex,"spindataset");
        oAttrib.put(cur);
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
        if(h5path != NULL) sprintf(hroot,"%s",(const char*)h5path);
      }
      cur=cur->next;
    }

    map<string,BsplineBasisType*>::iterator git(myBasis.find("0"));
    if(git == myBasis.end())
    {
      abasis = new OGT;
      if(OpenEndGrid)
        abasis->setGrid(LowerBox[0],UpperBox[0], LowerBox[1],UpperBox[1],LowerBox[2],UpperBox[2],
            BoxGrid[0],BoxGrid[1],BoxGrid[2]);
      else//need to offset for closed end
        abasis->setGrid(LowerBox[0],UpperBox[0],LowerBox[1],UpperBox[1],LowerBox[2],UpperBox[2],
            BoxGrid[0]-1,BoxGrid[1]-1,BoxGrid[2]-1);
      myBasis["0"]=abasis;
    } 
    else
    {
      abasis=dynamic_cast<OGT*>((*git).second);
    }

    abasis->setTwistAngle(targetPtcl.Lattice.k_cart(TwistAngle));
    typedef GroupedOrbitalSet<OGT>        SPOSetType;             
    StorageType inData(BoxGrid[0],BoxGrid[1],BoxGrid[2]);
    SPOSetType* psi= new SPOSetType(norb);
    bool complex2real = myParam->hasComplexData(hfileID);

#if defined(QMC_COMPLEX)
    if(!complex2real) 
    {
      app_error() << "  Real wavefunctions cannot be used with QMC_COMPLEX=1" << endl;
      abort(); //FIXABORT
    }
    complex2real = false;//reset to false
#endif

    if(complex2real)
    {
      Array<ComplexType,3> inTemp(BoxGrid[0],BoxGrid[1],BoxGrid[2]);
      for(int iorb=0; iorb<norb; iorb++) 
      {
        ostringstream wnshort;
        string eigvName=myParam->getEigVectorName(hroot,occSet[iorb]/degeneracy,spinIndex);
        wnshort<<curH5Fname << "#"<<occSet[iorb]/degeneracy << "#" << spinIndex;
        map<string,StorageType*>::iterator it(BigDataSet.find(wnshort.str()));
        if(it == BigDataSet.end()) {
          app_log() << "   Reading spline function " << eigvName << " (" << wnshort.str()  << ")" << endl;
          StorageType* newP =new StorageType;
          HDFAttribIO<Array<ComplexType,3> > dummy(inTemp);
          dummy.read(hfileID,eigvName.c_str());
          BLAS::copy(inTemp.size(),inTemp.data(),inData.data());
          BigDataSet[wnshort.str()]=newP;
          abasis->add(iorb,inData,newP);
        } else {
          app_log() << "   Reusing spline function " << eigvName << " (" << wnshort.str()  << ")" << endl;
          abasis->add(iorb,(*it).second);
        }
      } 
    }
    else 
    {
      for(int iorb=0; iorb<norb; iorb++) {
        ostringstream wnshort;
        string eigvName=myParam->getEigVectorName(hroot,occSet[iorb]/degeneracy,spinIndex);
        wnshort<<curH5Fname << "#"<<occSet[iorb]/degeneracy << "#" << spinIndex;
        map<string,StorageType*>::iterator it(BigDataSet.find(wnshort.str()));
        if(it == BigDataSet.end()) {
          app_log() << "   Reading spline function " << eigvName << " (" << wnshort.str()  << ")" << endl;
          StorageType* newP=new StorageType;
          HDFAttribIO<StorageType> dummy(inData);
          dummy.read(hfileID,eigvName.c_str());
          BigDataSet[wnshort.str()]=newP;
          abasis->add(iorb,inData,newP);
        } else {
          app_log() << "   Reusing spline function " << eigvName << " (" << wnshort.str()  << ")" << endl;
          abasis->add(iorb,(*it).second);
        }
      } 
    }
    psi->add(abasis);

    H5Fclose(hfileID);
    hfileID=-1;
    return psi;
  }

  /** create a SlaterDeterminant
   * @param cur xmlnode containing \<slaterdeterminant\>
   * @return a SlaterDeterminant
   *
   * @warning MultiSlaterDeterminant is not working yet.
   */
  SPOSetBase*
  TricubicBsplineSetBuilder::createSPOSet(xmlNodePtr cur){
    //return createSPOSetWithEG();
    //stage 1: get the general parameters, e.g., source, twist angle
    int norb(0);
    OhmmsAttributeSet aAttrib;
    aAttrib.add(norb,"orbitals"); aAttrib.add(norb,"size");
    aAttrib.add(curH5Fname,"href"); 
    aAttrib.put(cur);

    if(curH5Fname.empty())
    {
      app_error() << "No Valid HDF5 is provided for TricubicBsplineSetBuilder." << endl;
      app_error() << "Abort TricubicBsplineSetBuilder::createSPOSet(xmlNodePtr cur)" << endl;
      abort(); //FIXABORT
    }

    if(norb ==0) {
      app_error() << "TricubicBsplineSetBuilder::createSPOSet failed. Check the attribte orbitals." << endl;
      abort(); //FIXABORT
    }

    app_log() << "    HDF5 File = " << curH5Fname << endl;

    hfileID = H5Fopen(curH5Fname.c_str(),H5F_ACC_RDWR,H5P_DEFAULT);
    //check the version
    myParam->checkVersion(hfileID);
    //check htag multiple times
    myParam->put(rootNode);
    myParam->put(cur);

    string tname=myParam->getTwistAngleName();
    HDFAttribIO<PosType> hdfobj_twist(TwistAngle);
    hdfobj_twist.read(hfileID,tname.c_str());
    app_log() << "  Twist-angle " << TwistAngle << endl;

    RealType t2=dot(TwistAngle,TwistAngle);

    //here, we can also check if the cell is orthorohmbic or not 
    if(t2<numeric_limits<RealType>::epsilon()) //gamma point
    {
      app_log() << "  Gamma point. Using TricubicBsplineSet<ValueType> " << endl;
      TricubicBsplineSet<ValueType>* abasis=0;
      return createBsplineBasisSet(cur,abasis);
    }
    else
    {
#if defined(QMC_COMPLEX)
      app_log() << "  Non-gamma point. Using TricubicBsplineTwistSet<ValueType> " << endl;
      TricubicBsplineTwistSet<ValueType>* abasis=0;
      return createBsplineBasisSet(cur,abasis);
#else
      app_error() << "  Twist-angle other than Gamma point cannot work with real wavefunctions." << endl;
      abort(); //FIXABORT
      retrun 0;
#endif
    }
  }

  /** Create B-spline for the electron-gas wavefunctions
   *
   * This function is to test bspline against analytic electron-gas
   * single-particle orbitals.
   */
  SPOSetBase* TricubicBsplineSetBuilder::createSPOSetWithEG()
  {
    return 0; //disable
//#if defined(QMC_COMPLEX)
//    app_error() << " TricubicBsplineSetBuilder::createSPOSetWithEG cannot be used for QMC_COMPLEX=1" << endl;
//    OHMMS::Controller->abort();
//    return 0; //to make some compilers happy
//#else
//    int norb=7;
//    std::vector<int> npts(3);
//    npts[0]=BoxGrid[0]-1; 
//    npts[1]=BoxGrid[1]-1;
//    npts[2]=BoxGrid[2]-1;
//    StorageType inData(npts[0],npts[1],npts[2]);
//
//    RealType dx=(UpperBox[0]-LowerBox[0])/static_cast<RealType>(npts[0]);
//    RealType dy=(UpperBox[1]-LowerBox[1])/static_cast<RealType>(npts[1]);
//    RealType dz=(UpperBox[2]-LowerBox[2])/static_cast<RealType>(npts[2]);
//
//    SPOSetType* psi= new SPOSetType(norb);
//    map<string,OrbitalGroupType*>::iterator git(myBasis.find("0"));
//    OrbitalGroupType *abasis=0;
//    if(git == myBasis.end())
//    {
//      abasis = new OrbitalGroupType;
//      abasis->setGrid(LowerBox[0],UpperBox[0],LowerBox[1],UpperBox[1],LowerBox[2],UpperBox[2],
//          npts[0],npts[1],npts[2]);
//      myBasis["0"]=abasis;
//
//      int nc=1;
//      int nat=targetPtcl.getTotalNum();
//      int nup=nat/2;
//      HEGGrid<RealType,OHMMS_DIM> egGrid(targetPtcl.Lattice);
//      int nkpts=(nup-1)/2;
//
//      app_log() << "Number of kpoints " << nkpts << endl;
//      //create a E(lectron)G(as)O(rbital)Set
//      egGrid.createGrid(nc,nkpts);
//      RealEGOSet* eg=new RealEGOSet(egGrid.kpt,egGrid.mk2); 
//      char wfshortname[16];
//      for(int iorb=0; iorb<norb; iorb++) {
//        sprintf(wfshortname,"b%d",iorb);
//        for(int ix=0; ix<npts[0]; ix++) {
//          double x(dx*ix+LowerBox[0]);
//          for(int iy=0; iy<npts[1]; iy++) {
//            double y(dy*iy+LowerBox[1]);
//            double z(LowerBox[2]);
//            for(int iz=0; iz<npts[2]; iz++) {
//              inData(ix,iy,iz)=eg->f(PosType(x,y,z),iorb); z+=dz;
//            }
//          }
//        }
//        StorageType* newP=0;
//        map<string,StorageType*>::iterator it(BigDataSet.find(wfshortname));
//        if(it == BigDataSet.end()) {
//          newP=new StorageType;
//          BigDataSet[wfshortname]=newP;
//          abasis->add(iorb,inData,newP);
//          app_log() << "   Using spline function for EG " << wfshortname << endl;
//        } 
//      }
//      delete eg;
//    } 
//    else
//    {
//      abasis=(*git).second;
//    }
//
//    psi->add(abasis);
//    return psi;
//#endif
  }
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
