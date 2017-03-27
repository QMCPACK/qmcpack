//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
//This will go to cmake file
//template functions are slightly better with icpc/g++
//#define USE_BSPLINE_TEMP 1
#include "QMCWaveFunctions/PlaneWave/PWParameterSet.h"
#include "QMCWaveFunctions/Bspline3DSetTemp.h"
#include "QMCWaveFunctions/Bspline3DSet.h"
#include "QMCWaveFunctions/Bspline3DSetTrunc.h"
#include "QMCWaveFunctions/GroupedOrbitalSet.h"
#include "QMCWaveFunctions/TricubicBsplineSetBuilder.h"
#include "QMCWaveFunctions/OrbitalBuilderBase.h"
#include "OhmmsData/AttributeSet.h"
#include "ParticleIO/ParticleLayoutIO.h"
#include "Utilities/ProgressReportEngine.h"
#include "Message/CommOperators.h"
//#define DEBUG_BSPLINE_EG
namespace qmcplusplus
{

//initialize the static data member
//map<std::string,TricubicBsplineSetBuilder::StorageType*> TricubicBsplineSetBuilder::BigDataSet;
map<std::string,TricubicBsplineSetBuilder::RSOType*> TricubicBsplineSetBuilder::BigDataSet;

TricubicBsplineSetBuilder::TricubicBsplineSetBuilder(ParticleSet& p, PtclPoolType& psets, xmlNodePtr cur):
  targetPtcl(p),ptclPool(psets),OpenEndGrid(false),TranslateGrid(false),FloatingGrid(false),
  CurSPOSize(0),BoxGrid(2),BoxDup(1),LowerBox(0.0),UpperBox(1.0),
  rootNode(cur)
{
  ClassName="TricubicBsplineSetBuilder";
  basisLattice=p.Lattice;
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
bool TricubicBsplineSetBuilder::put(xmlNodePtr cur)
{
  ReportEngine PRE(ClassName,"put(xmlNodePtr)");
  //propagate the communicator
  myParam->initCommunicator(myComm);
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
    if(cname.find("cell")<cname.size())// == "supercell")
    {
      LatticeParser a(basisLattice);
      a.put(cur);
      for(int idim=0; idim<DIM; idim++)
        UpperBox[idim]=basisLattice.R(idim,idim);
    }
    else
      if(cname == "grid")
        //if(cname == "grid")
      {
        std::string closedEnd("yes");
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
        if(del>0 && del<UpperBox[idir])
          UpperBox[idir]=rf;
        LowerBox[idir]=ri;
        BoxGrid[idir]=npts;
        OpenEndGrid = (closedEnd != "yes");
      }
      else
        if(cname == "expand")//expand
        {
          putContent(BoxDup,cur);//
          myParam->BoxDup=BoxDup;
          TranslateGrid = true;
        }
    //else if(cname == "center")//move center of the grid
    //{
    //  std::string move("no");
    //  putContent(move,cur);
    //  TranslateGrid = (move =="yes");
    //}
    cur=cur->next;
  }
  PRE << "TricubicBsplineSetBuilder set a global box\n";
  PRE << "Lower box =" << LowerBox << "\n";
  PRE << "Upper box =" << UpperBox << "\n";
  PRE << "Box grid ="  << BoxGrid << "\n";
  PRE << "Box Expansion ="  << BoxDup << "\n";
  return true;
}

/** create a SingleParticleSetBase
 * @param cur xmlnode containing \<slaterdeterminant\>
 * @return a SlaterDeterminant
 *
 * @warning MultiSlaterDeterminant is not working yet.
 */
SPOSetBase*
TricubicBsplineSetBuilder::createSPOSet(xmlNodePtr cur)
{
  ReportEngine PRE(ClassName, "createSPOSet");
  //return createSPOSetWithEG();
  //stage 1: get the general parameters, e.g., source, twist angle
  int norb(0);
  std::string detname("0");
  OhmmsAttributeSet aAttrib;
  aAttrib.add(norb,"orbitals");
  aAttrib.add(norb,"size");
  aAttrib.add(curH5Fname,"href");
  aAttrib.add(detname,"id");
  aAttrib.put(cur);
  if(curH5Fname.empty())
  {
    //fatal
    PRE.error("No Valid HDF5 is provided for TricubicBsplineSetBuilder.",true);
  }
  if(norb ==0)
  {
    norb=targetPtcl.last(CurSPOSize)-targetPtcl.first(CurSPOSize);
    PRE <<"TricubicBsplineSetBuilder::createSPOSet missing size of the determinant.\n";
    PRE <<"Overwrite by the particle groups. " << norb << "\n";
    std::ostringstream norb_size;
    norb_size << norb;
    xmlNewProp(cur,(const xmlChar*)"size",(const xmlChar*)norb_size.str().c_str());
  }
  PRE <<"HDF5 File = " << curH5Fname << "\n";
  if(is_manager())
    hfileID = H5Fopen(curH5Fname.c_str(),H5F_ACC_RDWR,H5P_DEFAULT);
  else
    hfileID = -1;
  //check the version
  myParam->checkVersion(hfileID);
  //check htag multiple times
  myParam->put(rootNode);
  myParam->put(cur);
  if(is_manager())
  {
    std::string tname=myParam->getTwistAngleName();
    HDFAttribIO<PosType> hdfobj_twist(TwistAngle);
    hdfobj_twist.read(hfileID,tname.c_str());
  }
  myComm->bcast(TwistAngle);
  PRE<<"Twist-angle " << TwistAngle << "\n";
  //stage 2: analyze the supercell to choose a specialized bspline function
  bool atGamma= (dot(TwistAngle,TwistAngle)<std::numeric_limits<RealType>::epsilon());
  bool localize = myParam->Rcut>0.0;
  //RealType offdiag=0.0;
  //for(int idim=0; idim<DIM; idim++)
  //  for(int jdim=0; jdim<DIM; jdim++)
  //  {
  //    if(idim != jdim) offdiag+=std::abs(targetPtcl.Lattice.R(idim,jdim));
  //  }
  //bool orthorhombic=(offdiag< std::numeric_limits<RealType>::epsilon());
  bool orthorhombic=basisLattice.DiagonalOnly;
  //stage 3: check if the grid needs to be modified
  BoxDup=myParam->BoxDup;
  if(BoxDup[0]*BoxDup[1]*BoxDup[2]>1)
  {
    TranslateGrid=true;
  }
  else
  {
    TranslateGrid=(myParam->BufferRadius>0.0)? true:false;
  }
  //stage 4: create one of the derived classes from BsplineBasisType (Bspline3DSetBase)
  //constants for template instantiations
  const bool ORTHO=true;
  const bool NONORTHO=false;
  const bool TRUNC=true;
  const bool NOTRUNC=false;
  Timer t1;
  //determine activeBasis
  activeBasis=0;
  std::map<std::string,BsplineBasisType*>::iterator git(myBasis.find(detname));
  if(git == myBasis.end())
  {
    if(atGamma) //gamma point
    {
      if(localize)
      {
        if(TranslateGrid)
        {
          PRE <<"Bspline3DSet_MLW with truncated grid\n";
          activeBasis=new Bspline3DSet_MLW;
        }
        else
        {
          if(orthorhombic)
          {
            PRE <<"Bspline3DSet_Ortho_Trunc\n";
            activeBasis = new Bspline3DSet_Ortho_Trunc;
          }
          else
          {
            PRE <<"Bspline3DSet_Gen_Trunc\n";
            activeBasis = new Bspline3DSet_Gen_Trunc;
          }
        }
      }
      else
      {
        if(orthorhombic)
        {
          PRE <<"Bspline3DSet_Ortho\n";
          activeBasis = new Bspline3DSet_Ortho;
        }
        else
        {
          PRE <<"Bspline3DSet_Gen\n";
          activeBasis = new Bspline3DSet_Gen;
        }
      }
    }
#if defined(QMC_COMPLEX)
    else //any k-point
    {
      PRE <<"Bspline3DSet_Twist\n";
      activeBasis = new Bspline3DSet_Twist;
    }
#endif
    if(activeBasis == 0)
    {
      //fatal error
      PRE.error("Failed to create a Bspline functon.",true);
    }
    //initialize the grid
    initGrid();
    myBasis[detname]=activeBasis;
    activeBasis->setOrbitalSetSize(norb);
    //activeBasis->setLattice(targetPtcl.Lattice);
    activeBasis->setLattice(basisLattice);
    activeBasis->setTwistAngle(targetPtcl.Lattice.k_cart(TwistAngle));
    if(localize)
      activeBasis->setRcut(myParam->Rcut);
    //now parse everything else
    setBsplineBasisSet(cur);
  }
  else
  {
    //nothing to do, reuse an existing one
    activeBasis=(*git).second;
  }
  if(hfileID>-1)
  {
    H5Fclose(hfileID);
    hfileID=-1;
  }
  PRE << "Bspline Input = " << t1.elapsed() << " secs\n";
  //incremenet SPOSize
  CurSPOSize++;
  return activeBasis;
}

void TricubicBsplineSetBuilder::initGrid()
{
  if(is_manager())
  {
    //check the grid here
    herr_t status = H5Eset_auto(NULL, NULL);
    char gname[64];
    sprintf(gname,"/%s/grid",myParam->eigTag.c_str());
    status = H5Gget_objinfo (hfileID, gname, 0, NULL);
    if(status == 0)
    {
      hid_t g_in=H5Gopen(hfileID,gname);
      HDFAttribIO<TinyVector<int,DIM> > box_in(BoxGrid);
      box_in.read(g_in,"dimensions");
      PosType spacing;
      HDFAttribIO<PosType> spacing_in(spacing);
      spacing_in.read(g_in,"spacing");
      int yes=1;
      HDFAttribIO<int> tg(yes);
      tg.read(g_in,"translated");
      if(yes)
      {
        TranslateGrid=false;
        FloatingGrid=true;
      }
      tg.read(g_in,"closed");
      if(yes)
      {
        OpenEndGrid=false;
        UpperBox[0]=spacing[0]*(BoxGrid[0]-1);
        UpperBox[1]=spacing[1]*(BoxGrid[1]-1);
        UpperBox[2]=spacing[2]*(BoxGrid[2]-1);
      }
      else
      {
        OpenEndGrid=true;
        UpperBox[0]=spacing[0]*BoxGrid[0];
        UpperBox[1]=spacing[1]*BoxGrid[1];
        UpperBox[2]=spacing[2]*BoxGrid[2];
      }
    }
  }
  myComm->bcast(OpenEndGrid);
  myComm->bcast(BoxGrid);
  myComm->bcast(UpperBox);
  if(OpenEndGrid)
    activeBasis->setGrid(LowerBox[0],UpperBox[0], LowerBox[1],UpperBox[1],LowerBox[2],UpperBox[2],
                         BoxGrid[0],BoxGrid[1],BoxGrid[2],
                         basisLattice.BoxBConds[0],basisLattice.BoxBConds[1],basisLattice.BoxBConds[2],
                         OpenEndGrid);
  else//need to offset for closed end
    activeBasis->setGrid(LowerBox[0],UpperBox[0],LowerBox[1],UpperBox[1],LowerBox[2],UpperBox[2],
                         BoxGrid[0]-1,BoxGrid[1]-1,BoxGrid[2]-1,
                         basisLattice.BoxBConds[0],basisLattice.BoxBConds[1],basisLattice.BoxBConds[2],
                         OpenEndGrid);
  //copy the input grid
  dataKnot=activeBasis->bKnots;
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
//    app_error() << " TricubicBsplineSetBuilder::createSPOSetWithEG cannot be used for QMC_COMPLEX=1" << std::endl;
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
//    std::map<std::string,OrbitalGroupType*>::iterator git(myBasis.find("0"));
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
//      HEGGrid<RealType,DIM> egGrid(targetPtcl.Lattice);
//      int nkpts=(nup-1)/2;
//
//      app_log() << "Number of kpoints " << nkpts << std::endl;
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
//        std::map<std::string,StorageType*>::iterator it(BigDataSet.find(wfshortname));
//        if(it == BigDataSet.end()) {
//          newP=new StorageType;
//          BigDataSet[wfshortname]=newP;
//          abasis->add(iorb,inData,newP);
//          app_log() << "   Using spline function for EG " << wfshortname << std::endl;
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
