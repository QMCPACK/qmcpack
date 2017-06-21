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
    
    



#include "QMCApp/MO2Grid3D.h"
#include "QMCApp/ParticleSetPool.h"
#include "Utilities/OhmmsInfo.h"
#include "Particle/DistanceTableData.h"
#include "Particle/DistanceTable.h"
#include "ParticleIO/XMLParticleIO.h"
#include "QMCWaveFunctions/MolecularOrbitals/GridMolecularOrbitals.h"
#include "QMCWaveFunctions/LCOrbitals.h"
#include "OhmmsData/AttributeSet.h"
#include "Numerics/HDFSTLAttrib.h"
#include "Numerics/HDFTriCubicSpline.h"
#include "Utilities/Clock.h"
#include "QMCTools/QMCUtilities.h"

namespace qmcplusplus
{

MO2Grid3D::MO2Grid3D(int argc, char** argv): QMCAppBase(argc,argv),
  Electrons(0), Ions(0)
{
  ptclPool = new ParticleSetPool;
}

///destructor
MO2Grid3D::~MO2Grid3D()
{
  DEBUGMSG("MO2Grid3D::~MO2Grid3D")
  std::map<std::string,TriCubicSplineT<ValueType>* >::iterator it(SPOSet.begin());
  while(it != SPOSet.end())
  {
    if((*it).second)
      delete (*it).second;
    ++it;
  }
}

bool MO2Grid3D::validateXML()
{
  return true;
}

bool MO2Grid3D::execute()
{
  xmlNodePtr m_root = XmlDocStack.top()->getRoot();
  InFileRoot="generic";
  OhmmsAttributeSet pAttrib;
  pAttrib.add(InFileRoot,"id");
  pAttrib.add(InFileRoot,"name");
  pAttrib.put(m_root);
  //preserve the input order
  xmlNodePtr cur=m_root->children;
  xmlNodePtr wfsPtr=0;
  while(cur != NULL)
  {
    std::string cname((const char*)cur->name);
    if(cname == "particleset")
    {
      ptclPool->put(cur);
    }
    else
      if(cname == OrbitalBuilderBase::wfs_tag)
      {
        wfsPtr=cur;
      }
    cur=cur->next;
  }
  if(wfsPtr)
  {
    std::string ename("e");
    const xmlChar* t=xmlGetProp(wfsPtr,(const xmlChar*)"target");
    if(t)
    {
      ename=(const char*)t;
    }
    Electrons=ptclPool->getWalkerSet(ename);
    if(Electrons == 0)
    {
      ERRORMSG("ParticleSet with " << ename << " is not created. Abort")
      return false;
    }
    cur=wfsPtr->children;
    bool splitted=false;
    while(cur != NULL)
    {
      std::string cname((const char*)cur->name);
      if(cname == OrbitalBuilderBase::detset_tag)
      {
        dsetPtr=cur; //save determinantset
        std::string iname("i");
        const xmlChar* t=xmlGetProp(cur,(const xmlChar*)"source");
        if(t)
        {
          iname=(const char*)t;
        }
        Ions=ptclPool->getParticleSet(iname);
        if(Ions == 0)
        {
          ERRORMSG("ParticleSet with " << iname << " is not created. Abort")
          return false;
        }
        splitted = selectCore(cur);
        generateNumericalOrbitals(normalPtr);
      }
      cur=cur->next;
    }
    //remove the original determinantset
    xmlUnlinkNode(dsetPtr);
    xmlFreeNode(dsetPtr);
    //replace the original determinantset by ceorePtr
    xmlSetProp(corePtr,(const xmlChar*)"type",(const xmlChar*)"NumericalOrbital");
    xmlAddChild(wfsPtr,corePtr);
    std::string newfile=InFileRoot+".spline.xml";
    LOGMSG("New xml file " << newfile)
    XmlDocStack.top()->dump(newfile);
  }
  else
  {
    ERRORMSG("The input file does not contain wavefunction. Nothing to do.")
  }
  //return false to stop it
  return false;
}

/** temporary data to sort the atomicBasisSet and basisGroup
 */
struct BasisGroupType
{
  bool Excluded;
  int L;
  xmlNodePtr curPtr;
  inline BasisGroupType():Excluded(false),L(0),curPtr(0) { }
  inline BasisGroupType(xmlNodePtr cur):Excluded(false),L(0),curPtr(0)
  {
    put(cur);
  }
  inline bool put(xmlNodePtr cur)
  {
    curPtr=cur;
    std::string excluded("no");
    OhmmsAttributeSet pAttrib;
    pAttrib.add(excluded,"exclude");
    pAttrib.add(L,"l");
    pAttrib.put(cur);
    if(excluded == "yes")
      Excluded=true;
    return true;
  }
};

/** select basisGroup marked by excluded="yes"
 * @param cur determinantset node
 *
 * If there is any basisGroup marked by excluded="yes", two determinantset
 * nodes are created. Upon exit, corePtr contains new xml content
 * with type="NumericalOrbital".
 */
bool MO2Grid3D::selectCore(xmlNodePtr  cur)
{
  normalPtr=cur;
  typedef std::vector<BasisGroupType> AtomicBasisType;
  typedef std::map<int,AtomicBasisType* > RGroupType;
  RGroupType Rgroup;
  std::map<int,xmlNodePtr> RgroupPtr;
  SpeciesSet& ionSpecies(Ions->getSpeciesSet());
  xmlNodePtr basissetPtr=0;
  xmlNodePtr sdetPtr=0;
  xmlNodePtr splinePtr=0;
  cur=cur->children;
  while(cur != NULL)
  {
    std::string cname((const char*)cur->name);
    if(cname == OrbitalBuilderBase::basisset_tag)
    {
      basissetPtr=cur;//save the pointer
      xmlNodePtr cur1=cur->children;
      while(cur1 != NULL)
      {
        std::string cname1((const char*)cur1->name);
        if(cname1 == "atomicBasisSet")
        {
          int ionid
          =ionSpecies.addSpecies((const char*)xmlGetProp(cur1,(const xmlChar*)"elementType"));
          RGroupType::iterator rit(Rgroup.find(ionid));
          if(rit == Rgroup.end())
          {
            AtomicBasisType *acenter= new AtomicBasisType;
            xmlNodePtr cur2=cur1->children;
            while(cur2 != NULL)
            {
              std::string cname2((const char*)cur2->name);
              if(cname2 == "basisGroup")
              {
                acenter->push_back(BasisGroupType());
                acenter->back().put(cur2);
              }//basisGroup
              cur2=cur2->next;
            }
            Rgroup[ionid]=acenter;
            RgroupPtr[ionid]=cur1;
          }
        }//atomicBasisSet
        cur1=cur1->next;
      }//cur1
    }//basisset
    else
      if(cname == OrbitalBuilderBase::sd_tag)
      {
        sdetPtr=cur;//save slaterdeterminant
      }
      else
        if(cname == "cubicgrid")
        {
          splinePtr=cur;
        }
    cur=cur->next;
  }
  if(splinePtr == 0)
    //spline is missing, add one
  {
    std::vector<RealType> ri(3,-5.0);
    std::vector<RealType> rf(3,5.0);
    std::vector<int> npts(3,101);
    splinePtr = xmlNewNode(NULL,(const xmlChar*)"cubicgrid");
    for(int idir=0; idir<3; idir++)
    {
      xmlNodePtr x = xmlNewNode(NULL,(const xmlChar*)"grid");
      std::ostringstream dir_str, ri_str, rf_str, npts_str;
      dir_str << idir;
      ri_str << ri[idir];
      rf_str << rf[idir];
      npts_str << npts[idir];
      xmlNewProp(x,(const xmlChar*)"dir", (const xmlChar*)dir_str.str().c_str());
      xmlNewProp(x,(const xmlChar*)"ri",  (const xmlChar*)ri_str.str().c_str());
      xmlNewProp(x,(const xmlChar*)"rf",  (const xmlChar*)rf_str.str().c_str());
      xmlNewProp(x,(const xmlChar*)"npts",(const xmlChar*)npts_str.str().c_str());
      xmlAddChild(splinePtr,x);
    }
    xmlAddPrevSibling(basissetPtr,splinePtr);
  }
  std::vector<bool> mask;
  int offset=0, coreoffset=0;
  for(int iat=0; iat<Ions->getTotalNum(); iat++)
  {
    RGroupType::iterator rit(Rgroup.find(Ions->GroupID[iat]));
    AtomicBasisType &acenter(*((*rit).second));
    for(int i=0; i<acenter.size(); i++)
    {
      int l=acenter[i].L;
      bool excludeit=acenter[i].Excluded;
      for(int m=0; m<2*l+1; m++)
        mask.push_back(excludeit);
      if(excludeit)
        coreoffset+=2*l+1;
      offset += 2*l+1;
    }
  }
  //nothing to split, everything is on the grid
  if(coreoffset == 0)
  {
    LOGMSG("There is no localized basis set to exclude")
    corePtr = copyDeterminantSet(dsetPtr, splinePtr);
    return false;
  }
  //create two determinantset.
  //normalPtr is parsed by generateNumericalOrbitals
  //corePtr is used to write a new input file
  normalPtr = xmlCopyNode(dsetPtr,2);
  corePtr = xmlCopyNode(dsetPtr,2);
  xmlNodePtr core_1= xmlAddChild(corePtr,xmlCopyNode(splinePtr,1));
  xmlNodePtr normal_1 = xmlAddChild(normalPtr,xmlCopyNode(splinePtr,1));
  xmlNodePtr core_1_1 = xmlCopyNode(basissetPtr,2);
  xmlNodePtr normal_1_1 = xmlCopyNode(basissetPtr,2);
  RGroupType::iterator rit(Rgroup.begin()),rit_end(Rgroup.end());
  while(rit != rit_end)
  {
    AtomicBasisType &acenter(*((*rit).second));
    xmlNodePtr core_2=xmlCopyNode(RgroupPtr[(*rit).first],2);
    xmlNodePtr normal_2=xmlCopyNode(RgroupPtr[(*rit).first],2);
    for(int i=0; i<acenter.size(); i++)
    {
      if(acenter[i].Excluded)
      {
        xmlAddChild(core_2,xmlCopyNode(acenter[i].curPtr,1));
      }
      else
      {
        xmlAddChild(normal_2,xmlCopyNode(acenter[i].curPtr,1));
      }
    }
    xmlAddChild(core_1_1,core_2);
    xmlAddChild(normal_1_1,normal_2);
    delete (*rit).second;
    ++rit;
  }
  core_1=xmlAddSibling(core_1,core_1_1);
  normal_1=xmlAddSibling(normal_1,normal_1_1);
  core_1 = xmlAddSibling(core_1,xmlCopyNode(sdetPtr,2));
  normal_1 = xmlAddSibling(normal_1,xmlCopyNode(sdetPtr,2));
  cur=sdetPtr->children;
  while(cur != NULL)
  {
    std::string cname((const char*)cur->name);
    if(cname == OrbitalBuilderBase::det_tag)
    {
      //copy determinant
      xmlNodePtr core_2=xmlCopyNode(cur,2);
      xmlNodePtr normal_2=xmlCopyNode(cur,2);
      int norb=atoi((const char*)xmlGetProp(cur,(const xmlChar*)"orbitals"));
      std::vector<RealType> tmpC(offset*offset), occ(offset), C(norb*offset);
      Matrix<RealType> A(norb,coreoffset), B(norb,offset-coreoffset);
      bool foundCoeff=false;
      xmlNodePtr cur1=cur->children;
      while(cur1!=NULL)
      {
        std::string cname1((const char*)cur1->name);
        if(cname1 == "occupation")
        {
          putContent(occ.begin(), occ.end(),cur1);
        }
        else
          if(cname1 == "coefficient" || cname1 == "parameter")
          {
            foundCoeff = true;
            putContent(tmpC.begin(), tmpC.end(),cur1);
          }
        cur1=cur1->next;
      }
      if(foundCoeff)
      {
        int iorb(0),n(0);
        std::vector<ValueType>::iterator tit(tmpC.begin());
        std::vector<ValueType>::iterator cit(C.begin());
        while(iorb<norb)
        {
          if(occ[n]>std::numeric_limits<RealType>::epsilon())
          {
            copy(tit,tit+offset,cit);
            iorb++;
            cit+=offset;
          }
          ++n;
          tit+=offset;
        }
        int kk=0;
        for(iorb=0; iorb<norb; iorb++)
        {
          int ka(0), kb(0);
          for(int k=0; k<offset; k++, kk++)
          {
            if(mask[k])
              A(iorb,ka++)=C[kk];
            else
              B(iorb,kb++)=C[kk];
          }
        }
        getEigVectors(core_2,A);
        getEigVectors(normal_2,B);
      }
      //add basis
      xmlAddChild(core_2,copyDeterminant(cur,true));
      core_2 = xmlAddChild(core_1,core_2);
      normal_2 = xmlAddChild(normal_1,normal_2);
    }//determinant
    cur=cur->next;
  }
  return true;
}

/** copy determinant node for numerical orbitals
 * @param cur xmlnode to be copied
 * @param addg if ture, append _g to the id attribute
 * @return a new xmlnode for printout
 *
 * A new node \<spline\> is created to hold the information about the numerical orbitals.
 * The attribute basis="cubicgrid" is added to point to the cubicgrid this
 * spline node will use.
 * The root name of hdf5 files for each orbital is set to src attribute.
 */
xmlNodePtr MO2Grid3D::copyDeterminant(xmlNodePtr cur, bool addg)
{
  xmlNodePtr aptr_hdf5 = xmlNewNode(NULL,(const xmlChar*)"spline");
  xmlNewProp(aptr_hdf5,(const xmlChar*)"basis",(const xmlChar*)"cubicgrid");
  char oname[128];
  const xmlChar* refptr=xmlGetProp(cur,(const xmlChar*)"ref");
  if(refptr)
  {
    sprintf(oname,"%s.%s",InFileRoot.c_str(),(const char*)refptr);
  }
  else
  {
    refptr=xmlGetProp(cur,(const xmlChar*)"id");
    sprintf(oname,"%s.%s",InFileRoot.c_str(),(const char*)refptr);
  }
  xmlNewProp(aptr_hdf5,(const xmlChar*)"src",(const xmlChar*)oname);
  return aptr_hdf5;
}

/** copy determinantset node
 * @param cur original determinantset node
 * @param splinePtr xml node containing cubicgrid
 * @return new determinantset node for printout
 *
 * This builds a xml node for the new input file containing numerical orbitals.
 */
xmlNodePtr MO2Grid3D::copyDeterminantSet(xmlNodePtr cur, xmlNodePtr splinePtr)
{
  xmlNodePtr acopy = xmlCopyNode(cur,2);
  xmlNodePtr next1 = xmlAddChild(acopy,xmlCopyNode(splinePtr,1));
  cur=cur->children;
  while(cur != NULL)
  {
    std::string cname((const char*)cur->name);
    if(cname == OrbitalBuilderBase::sd_tag)
    {
      next1 = xmlAddSibling(next1,xmlCopyNode(cur,2));
      xmlNodePtr cur1=cur->children;
      while(cur1 != NULL)
      {
        std::string cname1((const char*)cur1->name);
        if(cname1 == OrbitalBuilderBase::det_tag)
        {
          xmlNodePtr newdet=xmlCopyNode(cur1,2);
          xmlAddChild(newdet,copyDeterminant(cur1,false));
          xmlAddChild(next1,newdet);
        }
        cur1=cur1->next;
      }
    }
    cur=cur->next;
  }
  return acopy;
}

void
MO2Grid3D::getEigVectors(xmlNodePtr adet, const Matrix<RealType>& A)
{
  std::ostringstream eig,b_size;
  b_size << A.rows();
  eig.setf(std::ios::scientific, std::ios::floatfield);
  eig.setf(std::ios::right,std::ios::adjustfield);
  eig.precision(14);
  eig << "\n";
  int btot=A.size(),b=0;
  int n=btot/4;
  int dn=btot-n*4;
  for(int k=0; k<n; k++)
  {
    eig << std::setw(22) << A(b) << std::setw(22) << A(b+1)
        << std::setw(22) << A(b+2) << std::setw(22) <<  A(b+3) << "\n";
    b += 4;
  }
  for(int k=0; k<dn; k++)
  {
    eig << std::setw(22) << A(b);
  }
  xmlNodePtr det_data
  = xmlNewTextChild(adet,NULL,(const xmlChar*)"coefficient",(const xmlChar*)eig.str().c_str());
  xmlNewProp(det_data,(const xmlChar*)"size",(const xmlChar*)b_size.str().c_str());
  xmlNewProp(det_data,(const xmlChar*)"basisset",(const xmlChar*)"mo");
}

void
MO2Grid3D::copyOrbitalSet(std::map<std::string,TriCubicSplineT<ValueType>* >& other)
{
  std::map<std::string,TriCubicSplineT<ValueType>* >::iterator it(SPOSet.begin());
  while(it != SPOSet.end())
  {
    other[(*it).first]=(*it).second;
    ++it;
  }
  //clear the set so that it does not delete the objects
  SPOSet.clear();
}

xmlNodePtr MO2Grid3D::generateNumericalOrbitals(xmlNodePtr cur)
{
  LOGMSG("MO2Grid3D::generateNumericalOrbitals from Molecular Orbitals")
  TrialWaveFunction Psi;
  GridMolecularOrbitals Original((*Electrons),Psi,(*Ions));
  typedef GridMolecularOrbitals::BasisSetType BasisSetType;
  typedef LCOrbitals<BasisSetType>            SPOSetType;
  int nels=0;
  std::vector<RealType> ri(3,-5.0);
  std::vector<RealType> rf(3,5.0);
  std::vector<int> npts(3,101);
  std::vector<SPOSetType*> InOrbs;
  std::map<std::string,int> DetCounter;
  std::vector<xmlNodePtr> SlaterDetPtr;
  xmlNodePtr splinePtr=0;
  xmlNodePtr firstSlaterDetPtr=0;
  BasisSetType *basisSet=0;
  const xmlChar* a;
  std::vector<xmlNodePtr> nodeToBeRemoved;
  xmlNodePtr curRoot=cur;
  cur = cur->xmlChildrenNode;
  int idir=0;
  //first pass to check if basisset is used
  while(cur != NULL)
  {
    std::string cname((const char*)(cur->name));
    if(cname == "cubicgrid")
    {
      splinePtr=cur; // save spline data
      xmlNodePtr tcur = cur->xmlChildrenNode;
      while(tcur != NULL)
      {
        std::string tname((const char*)(tcur->name));
        if(tname == "grid")
        {
          a=xmlGetProp(tcur,(const xmlChar*)"dir");
          if(a)
          {
            idir=atoi((const char*)a);
          }
          a=xmlGetProp(tcur,(const xmlChar*)"ri");
          if(a)
          {
            ri[idir]=atof((const char*)a);
          }
          a=xmlGetProp(tcur,(const xmlChar*)"rf");
          if(a)
          {
            rf[idir]=atof((const char*)a);
          }
          a=xmlGetProp(tcur,(const xmlChar*)"npts");
          if(a)
          {
            npts[idir]=atoi((const char*)a);
          }
        }
        tcur = tcur->next;
      }
    }
    else
      if(cname == OrbitalBuilderBase::basisset_tag)
      {
        //basisptr=cur;
        nodeToBeRemoved.push_back(cur);
        basisSet = Original.addBasisSet(cur);
      }
    cur = cur->next;
  }
  //basisset is not given.  return 0
  if(basisSet == 0)
    return 0;
  cur = curRoot->xmlChildrenNode;
  while(cur != NULL)
  {
    std::string cname((const char*)(cur->name));
    if(cname == OrbitalBuilderBase::sd_tag)
    {
      if(firstSlaterDetPtr==0)
        firstSlaterDetPtr=cur;
      SlaterDetPtr.push_back(cur);
      nels=0;
      xmlNodePtr tcur = cur->xmlChildrenNode;
      while(tcur != NULL)
      {
        std::string tname((const char*)(tcur->name));
        if(tname == OrbitalBuilderBase::det_tag)
        {
          bool newset=true;
          std::string detname("det");
          a=xmlGetProp(tcur,(const xmlChar*)"id");
          if(a)
          {
            detname = (const char*)a;
          }
          else
          {
            //when id is not given, newset is set to false
            if(DetCounter.size())
              newset=false;
          }
          std::string detref(detname);
          a=xmlGetProp(tcur,(const xmlChar*)"ref");
          if(a)
          {
            detref=(const char*)a;
            if(DetCounter.find(detref) != DetCounter.end())
            {
              newset=false;
            }
          }
          a=xmlGetProp(tcur,(const xmlChar*)"orbitals");
          int norbs = atoi((const char*)a);
          xmlNodePtr c=tcur->children;
          while(c != NULL)
          {
            nodeToBeRemoved.push_back(c);
            c=c->next;
          }
          if(newset)
          {
            SPOSetType* psi=new SPOSetType(basisSet,nels);
            psi->put(tcur);
            psi->setName(detname);
            psi->resize(norbs);
            DetCounter[detname]=InOrbs.size();
            InOrbs.push_back(psi);
            //InOrbs[detname]=psi;
          }
          //add src attribute
          std::string detsrc=InFileRoot+"."+detref;
          if(xmlHasProp(tcur,(const xmlChar*)"src"))
          {
            xmlSetProp(tcur,(const xmlChar*)"src",(const xmlChar*)detsrc.c_str());
          }
          else
          {
            xmlNewProp(tcur,(const xmlChar*)"src",(const xmlChar*)detsrc.c_str());
          }
          nels+=norbs;
        }
        tcur = tcur->next;
      }
    }
    cur = cur->next;
  }
  //resize with respect to the number of electrons
  basisSet->resize(nels);
  //Need only one electron to calculate this
  //map<std::string,SPOSetType*>::iterator oit(InOrbs.begin());
  //Create one-dimensional grids for three orthogonal directions
  typedef LinearGrid<double> GridType;
  GridType *gridX=new GridType;
  GridType *gridY=new GridType;
  GridType *gridZ=new GridType;
  gridX->set(ri[0],rf[0],npts[0]);
  gridY->set(ri[1],rf[1],npts[1]);
  gridZ->set(ri[2],rf[2],npts[2]);
  typedef TriCubicSplineT<ValueType> NOType;
  XYZCubicGrid<RealType> *grid3 = new XYZCubicGrid<RealType>(gridX,gridY,gridZ);
  int ntot = npts[0]*npts[1]*npts[2];
  //vector<ValueType> phi(inorb->numOrbitals(),0.0);
  std::vector<ValueType> dat(ntot,0.0);
  Pooma::Clock timer;
  char oname[128];
  std::cout << "XYZCubicGrid " << std::endl;
  std::cout << " x " << ri[0] << " " << rf[0] << " " << npts[0] << std::endl;
  std::cout << " y " << ri[1] << " " << rf[1] << " " << npts[1] << std::endl;
  std::cout << " z " << ri[2] << " " << rf[2] << " " << npts[2] << std::endl;
  std::vector<RealType> lapsed_time(5,0.0);
  PosType pos(Electrons->R[0]);
  std::map<std::string,int>::iterator dit(DetCounter.begin());
  //loop over unique determinant sets
  while(dit != DetCounter.end())
  {
    SPOSetType* inorb=InOrbs[(*dit).second];
    std::string detID((*dit).first);
    for(int iorb=0; iorb<inorb->numOrbitals(); iorb++)
      //evaluate the values on the grid points
    {
      for(int ix=1; ix<npts[0]-1; ix++)
      {
        RealType x((*gridX)(ix));
        for(int iy=1; iy<npts[1]-1; iy++)
        {
          RealType y((*gridY)(iy));
          int ixyz=1+npts[2]*(iy+npts[1]*ix);
          for(int iz=1; iz<npts[2]-1; iz++, ixyz++)
          {
            PosType dr(x,y,(*gridZ)(iz));
            PosType newpos(Electrons->makeMove(0,dr-pos));
            dat[ixyz] = inorb->evaluate(*Electrons,0,iorb);
          }
        }
      }
      sprintf(oname,"%s.%s.wf%04d.h5",InFileRoot.c_str(),detID.c_str(),iorb);
      hid_t h_file = H5Fcreate(oname,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
      HDFAttribIO<std::vector<double> > dump(dat,npts);
      dump.write(h_file,"Orbital");
      H5Fclose(h_file);
    }
//      for(int iorb=0; iorb<inorb->numOrbitals(); iorb++) { //evaluate the values on the grid points
//        NOType *torb = new NOType(grid3);
//        timer.start();
//        for(int ix=1; ix<npts[0]-1; ix++) {
//          RealType x((*gridX)(ix));
//          for(int iy=1; iy<npts[1]-1; iy++) {
//            RealType y((*gridY)(iy));
//            for(int iz=1; iz<npts[2]-1; iz++) {
//              PosType dr(x,y,(*gridZ)(iz));
//              PosType newpos(Electrons->makeMove(0,dr-pos));
//              (*torb)(ix,iy,iz) = inorb->evaluate(*Electrons,0,iorb);
//            }
//          }
//        }
//        timer.stop();
//        lapsed_time[0]+= timer.cpu_time();
//
//        //spline
//        timer.start();
//        torb->reset(false);
//        timer.stop();
//        lapsed_time[1]+= timer.cpu_time();
//
//        //test
//        timer.start();
//        for(int ix=0; ix<npts[0]-1; ix++) {
//          double x((*gridX)(ix));
//          for(int iy=0; iy<npts[1]-1; iy++) {
//            double y((*gridY)(iy));
//            int offset=npts[2]*(iy+npts[1]*ix);
//            for(int iz=0; iz<npts[2]-1; iz++,offset++) {
//               TinyVector<double,3> p(x,y,(*gridZ)(iz));
//               dat[offset]=torb->evaluate(p);
//            }
//          }
//        }
//        timer.stop();
//        lapsed_time[2]+= timer.cpu_time();
//
//        timer.start();
//        sprintf(oname,"%s.%s.wf%04d.h5",InFileRoot.c_str(),detID.c_str(),iorb);
//        hid_t h_file = H5Fcreate(oname,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
//        HDFAttribIO<std::vector<double> > dump(dat,npts);
//        dump.write(h_file,"Orbital");
//        //HDFAttribIO<TriCubicSplineT<double> > dump1(*torb);
//        //dump1.write(h_file,"CubicSpline");
//        H5Fclose(h_file);
//        timer.stop();
//        lapsed_time[3]+= timer.cpu_time();
//
//        //add to SPOSet
//        sprintf(oname,"%s.%s.wf%04d",InFileRoot.c_str(),detID.c_str(),iorb);
//        SPOSet[oname]=torb;
//      }
    ++dit;
  }
  std::cout << "Timing results in sec" << std::endl;
  std::cout << "Function evaluation " << nels << " orbitals = " << lapsed_time[0] << std::endl;
  std::cout << "Spline coefficients = " << lapsed_time[1] << std::endl;
  std::cout << "Testing spline      = " << lapsed_time[2] << std::endl;
  std::cout << "Writing hdf5 files  = " << lapsed_time[3] << std::endl;
  //clean up temporary orbitals on the radial grid
  for(int i=0; i<InOrbs.size(); i++)
    delete InOrbs[i];
  //return the modified node
  return curRoot;
}
}
