//////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
//This will go to cmake file
//template functions are slightly better with icpc/g++
//#define USE_BSPLINE_TEMP 1
#include "QMCWaveFunctions/PlaneWave/PWParameterSet.h"
#include "QMCWaveFunctions/Bspline3DSetTemp.h"
#include "QMCWaveFunctions/Bspline3DSet.h"
#include "QMCWaveFunctions/GroupedOrbitalSet.h"
#include "QMCWaveFunctions/TricubicBsplineSetBuilder.h"
#include "QMCWaveFunctions/OrbitalBuilderBase.h"
#include "OhmmsData/AttributeSet.h"
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
      print_log = OrbitalBuilderBase::print_level>0;
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
    if(print_log)
    {
      app_log() << "  TricubicBsplineSetBuilder set a global box " << endl;
      app_log() << "    Lower box " << LowerBox << endl;
      app_log() << "    Upper box " << UpperBox << endl;
      app_log() << "    Box grid "  << BoxGrid << endl;
    }
    return true;
  }

  /** create a SingleParticleSetBase
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

    if(print_log) app_log() << "    HDF5 File = " << curH5Fname << endl;
    hfileID = H5Fopen(curH5Fname.c_str(),H5F_ACC_RDWR,H5P_DEFAULT);
    //check the version
    myParam->checkVersion(hfileID);
    //check htag multiple times
    myParam->put(rootNode);
    myParam->put(cur);

    string tname=myParam->getTwistAngleName();
    HDFAttribIO<PosType> hdfobj_twist(TwistAngle);
    hdfobj_twist.read(hfileID,tname.c_str());
    if(print_log)  app_log() << "    Twist-angle " << TwistAngle << endl;

    //stage 2: analyze the supercell to choose a specialized bspline function
    bool atGamma= (dot(TwistAngle,TwistAngle)<numeric_limits<RealType>::epsilon());
    bool localize = myParam->Rcut>0.0;
    RealType offdiag=0.0;
    for(int idim=0; idim<OHMMS_DIM; idim++)
      for(int jdim=0; jdim<OHMMS_DIM; jdim++)
      {
        if(idim != jdim) offdiag+=abs(targetPtcl.Lattice.R(idim,jdim));
      }
    bool orthorhombic=(offdiag< numeric_limits<RealType>::epsilon());

    //stage 3: create one of the derived classes from BsplineBasisType (Bspline3DSetBase)
    //constants for template instantiations
    const bool ORTHO=true;
    const bool NONORTHO=false;
    const bool TRUNC=true;
    const bool NOTRUNC=false;

    //determine activeBasis
    activeBasis=0;
    map<string,BsplineBasisType*>::iterator git(myBasis.find("0"));
    if(git == myBasis.end())
    {
      if(atGamma) //gamma point
      {
        if(localize)
        {
          if(orthorhombic)
          {
#if defined(USE_BSPLINE_TEMP)
            if(print_log)  app_log() << "    Bspline3DSet<ORTHO=true,TRUC=true> " << endl;
            activeBasis = new Bspline3DSet<ORTHO,TRUNC>;
#else
            if(print_log)  app_log() << "    Bspline3DSet_Ortho_Trunc" << endl;
            activeBasis = new Bspline3DSet_Ortho_Trunc;
#endif
          }
          else
          {
#if defined(USE_BSPLINE_TEMP)
            if(print_log) app_log() << "    Bspline3DSet<ORTHO=false,TRUC=true> " << endl;
            activeBasis = new  Bspline3DSet<NONORTHO,TRUNC>;
#else
            if(print_log) app_log() << "    Bspline3DSet_Gen_Trunc " << endl;
            activeBasis = new Bspline3DSet_Gen_Trunc;
#endif
          } 
        }
        else
        {
          if(orthorhombic)
          {
#if defined(USE_BSPLINE_TEMP)
            if(print_log)  app_log() << "    Bspline3DSet<ORTHO=true,TRUC=false> " << endl;
            activeBasis = new Bspline3DSet<ORTHO,NOTRUNC>;
#else
            if(print_log)  app_log() << "    Bspline3DSet_Ortho " << endl;
            activeBasis = new Bspline3DSet_Ortho;
#endif
          }
          else
          {
#if defined(USE_BSPLINE_TEMP)
            if(print_log)  app_log() << "    Bspline3DSet<ORTHO=false,TRUC=false> " << endl;
            activeBasis = new Bspline3DSet<NONORTHO,NOTRUNC>;
#else
            if(print_log)  app_log() << "    Bspline3DSet_Gen " << endl;
            activeBasis = new Bspline3DSet_Gen;
#endif
          }
        }
      }
#if defined(QMC_COMPLEX)
      else //any k-point
      {
        if(print_log) app_log() << "    Bspline3DSet_Twist" << endl;
        activeBasis = new Bspline3DSet_Twist;
      }
#endif

      if(activeBasis == 0) 
      {
        app_error() << "  Failed to create a Bspline functon. Abort. " << endl;
        abort(); //FIXABORT
      }

      if(OpenEndGrid)
        activeBasis->setGrid(LowerBox[0],UpperBox[0], LowerBox[1],UpperBox[1],LowerBox[2],UpperBox[2],
            BoxGrid[0],BoxGrid[1],BoxGrid[2]);
      else//need to offset for closed end
        activeBasis->setGrid(LowerBox[0],UpperBox[0],LowerBox[1],UpperBox[1],LowerBox[2],UpperBox[2],
            BoxGrid[0]-1,BoxGrid[1]-1,BoxGrid[2]-1);

      myBasis["0"]=activeBasis;
      activeBasis->setOrbitalSetSize(norb);
      activeBasis->setLattice(targetPtcl.Lattice);
      activeBasis->setTwistAngle(targetPtcl.Lattice.k_cart(TwistAngle));
      if(localize) activeBasis->setRcut(myParam->Rcut);
    } 
    else
    {
      activeBasis=(*git).second;
    }

    //now parse everything else
    setBsplineBasisSet(cur);

    H5Fclose(hfileID);
    hfileID=-1;

    return activeBasis;
  }

//  /** Create B-spline for the electron-gas wavefunctions
//   *
//   * This function is to test bspline against analytic electron-gas
//   * single-particle orbitals.
//   */
//  SPOSetBase* TricubicBsplineSetBuilder::createSPOSetWithEG()
//  {
//    return 0; //disable
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
//  }
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
