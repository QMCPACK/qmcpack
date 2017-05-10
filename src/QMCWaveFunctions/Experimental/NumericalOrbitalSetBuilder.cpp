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
    
    
#include "QMCWaveFunctions/NumericalOrbitalSetBuilder.h"
#include "OhmmsData/AttributeSet.h"
#include "Numerics/HDFSTLAttrib.h"
//#include "Numerics/HDFTriCubicSpline.h"

namespace qmcplusplus
{

NumericalOrbitalSetBuilder::NumericalOrbitalSetBuilder(ParticleSet& p,
    TrialWaveFunction& psi, PtclPoolType& psets):
  OrbitalBuilderBase(p,psi),
  LocalizedBasisBuilder(0),
  ptclPool(psets),
  GridXYZ(0), MOBasisSet(0)
{
}

/** add an antisymmetric many-body function to a trial wave function
 */
bool NumericalOrbitalSetBuilder::put(xmlNodePtr cur)
{
  std::string source("i");
  OhmmsAttributeSet aAttrib;
  aAttrib.add(source,"source");
  aAttrib.put(cur);
  //ionic system for the molecualr orbitals
  ParticleSet* ions=0;
  //initialize with the source tag
  PtclPoolType::iterator pit(ptclPool.find(source));
  if(pit != ptclPool.end())
  {
    LOGMSG("Molecular orbital uses an ionic system " << source)
    ions=(*pit).second;
  }
  if(LocalizedBasisBuilder == 0)
  {
    LocalizedBasisBuilder = new GridMolecularOrbitals(targetPtcl,
        targetPsi,*ions);
  }
  xmlNodePtr gridPtr=0;
  std::vector<xmlNodePtr> sdetPtr;
  cur=cur->children;
  while(cur != NULL)
  {
    std::string cname((const char*)(cur->name));
    if(cname == "cubicgrid")
    {
      gridPtr=cur;
    }
    else
      if(cname == sd_tag)
      {
        sdetPtr.push_back(cur);
      }
      else
        if(cname == basisset_tag)
        {
          MOBasisSet = LocalizedBasisBuilder->addBasisSet(cur);
        }
    cur = cur->next;
  }
  if(gridPtr == 0)
  {
    ERRORMSG("NumericalOrbitalSetBuilder::put failed. cubicgrid is not defined.")
    return false;
  }
  if(GridXYZ == 0)
  {
    GridXYZ = new XYZCubicGrid<RealType>;
    GridXYZ->put(gridPtr);
  }
  //For now, only one slaterdeterminant is used
  //add slaterdeterminant to a trial wave function
  SlaterDeterminant_t* sdet = addSlaterDeterminant(sdetPtr[0]);
  SlaterDetSet.push_back(sdet);
  targetPsi.addOrbital(sdet);
  //BasisSet needs to be resized
  if(MOBasisSet)
    MOBasisSet->resize(targetPtcl.getTotalNum());
  XMLReport("Creating a SlaterDeterminant wavefunction")
  return true;
}

/** create a SlaterDeterminant
 * @param cur xmlnode containing \<slaterdeterminant\>
 * @return a SlaterDeterminant
 *
 * @warning MultiSlaterDeterminant is not working yet.
 */
NumericalOrbitalSetBuilder::SlaterDeterminant_t*
NumericalOrbitalSetBuilder::addSlaterDeterminant(xmlNodePtr cur)
{
  std::string nogood("invalid");
  int is=0;
  SlaterDeterminant_t* sdet=new SlaterDeterminant_t;
  int first = 0;
  cur = cur->xmlChildrenNode;
  while(cur != NULL)
  {
    std::string cname((const char*)(cur->name));
    if(cname == det_tag)
    {
      std::string id("det"),ref(nogood);
      int norb(0);
      OhmmsAttributeSet aAttrib;
      aAttrib.add(id,"id");
      aAttrib.add(ref,"ref");
      aAttrib.add(norb,"orbitals");
      aAttrib.put(cur);
      //when ref is not given, overwrite ref by id
      if(ref == nogood)
        ref=id;
      SPOSetType* psi= 0;
      //map<std::string,SPOSetType*>::iterator sit(SPOSet.find(id));
      std::map<std::string,SPOSetType*>::iterator sit(SPOSet.find(ref));
      if(sit == SPOSet.end())
      {
        LOGMSG("Adding a single-particle orbital " << id)
        psi = createSPOSet(cur,ref,norb);
      }
      else
      {
        LOGMSG("Reusing an existing single-particle orbital " << id)
        psi = (*sit).second;
      }
      std::map<std::string,Det_t*>::iterator dit(DetSet.find(id));
      //check if determinant's id is valid
      if(dit != DetSet.end() && SlaterDetSet.empty())
      {
        std::ostringstream idassigned(id);
        idassigned << is;
        ERRORMSG("determinants cannot have same id. id will be assigned to " << id)
      }
      Det_t *adet =0;
      dit = DetSet.find(id);
      if(dit == DetSet.end())
        //need to add a new Determinant
      {
        std::cout << "Checking determinant , first, norb " << std::endl;
        adet = new Det_t(*psi,first);
        adet->set(first,norb);
        first+=norb;
        DetSet[id]=adet;
        LOGMSG("Adding a determinant " << id)
      }
      else
      {
        adet=(*dit).second;
      }
      //add determinant to slaterdeterminant
      sdet->add(adet);
    }
    cur = cur->next;
  }
  return sdet;
}

/** create a SingleParticleOrbitalSet containng norb orbitals
 * @param cur xmlnode to process
 * @param ref name of SPO set
 * @param norb number of single-particle orbitals
 * @return a new MixedSPOSet
 */
NumericalOrbitalSetBuilder::SPOSetType*
NumericalOrbitalSetBuilder::createSPOSet(xmlNodePtr cur, const std::string& ref, int norb)
{
  std::vector<int> npts(3);
  npts[0]=GridXYZ->nX;
  npts[1]=GridXYZ->nY;
  npts[2]=GridXYZ->nZ;
  std::vector<RealType> inData(npts[0]*npts[1]*npts[2]);
  int first=LocalizedOrbitals.size();
  SPOSetType* psi= new SPOSetType;
  cur=cur->children;
  while(cur != NULL)
  {
    std::string cname((const char*)(cur->name));
    if(cname == "coefficient" || cname == "parameter")
    {
      LOType *phi=0;
      std::map<std::string,LOType*>::iterator it(LocalizedOrbitals.find(ref));
      if(it == LocalizedOrbitals.end())
      {
        phi=new LOType(MOBasisSet,first);
        phi->put(cur->parent);
        phi->setName(ref);
        LocalizedOrbitals[ref]=phi;
      }
      else
      {
        phi = (*it).second;
      }
      psi->setLocalizedOrbitals(phi);
    }
    else
      if(cname == "spline")
      {
        const char* hroot = (const char*)xmlGetProp(cur,(const xmlChar*)"src");
        char wfname[128], hfile[128];
        for(int iorb=0; iorb<norb; iorb++)
        {
          sprintf(wfname,"%s.wf%04d",hroot,iorb);
          sprintf(hfile,"%s.wf%04d.h5",hroot,iorb);
          NumericalOrbitalType* neworb=0;
          std::map<std::string,NumericalOrbitalType*>::iterator it(NumericalOrbitals.find(wfname));
          if(it == NumericalOrbitals.end())
          {
            neworb=new NumericalOrbitalType(GridXYZ);
            hid_t h_file = H5Fopen(hfile,H5F_ACC_RDWR,H5P_DEFAULT);
            HDFAttribIO<std::vector<RealType> > dummy(inData,npts);
            dummy.read(h_file,"/Orbital");
            H5Fclose(h_file);
            //////////////////////////////
            //HDFAttribIO<NumericalOrbitalType> dummy(*neworb);
            //dummy.read(h_file,"CubicSpline");
            ////////////////////////////
            neworb->reset(inData.begin(), inData.end(), false);
            NumericalOrbitals[wfname]=neworb;
          }
          else
          {
            neworb = (*it).second;
          }
          psi->add(neworb);
        }
      }
    cur=cur->next;
  }
  return psi;
}
}
