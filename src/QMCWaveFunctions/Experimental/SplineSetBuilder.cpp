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
    
    
#include "OhmmsData/AttributeSet.h"
#include "Numerics/HDFSTLAttrib.h"
#include "QMCWaveFunctions/SplineSetBuilder.h"
//#include "Numerics/HDFTriCubicSpline.h"

namespace qmcplusplus
{

SplineSetBuilder::SplineSetBuilder(ParticleSet& p, PtclPoolType& psets):
  targetPtcl(p),ptclPool(psets),GridXYZ(0)
{
}

/** Initialize the cubic grid
 */
bool SplineSetBuilder::put(xmlNodePtr cur)
{
  if(GridXYZ==0)
  {
    GridXYZ = new GridType;
    GridXYZ->put(cur);
  }
  return true;
}

/** create a SlaterDeterminant
 * @param cur xmlnode containing \<slaterdeterminant\>
 * @return a SlaterDeterminant
 *
 * @warning MultiSlaterDeterminant is not working yet.
 */
SPOSetBase*
SplineSetBuilder::createSPOSet(xmlNodePtr cur)
{
  std::string hrefname("NONE");
  int norb(0);
  int degeneracy(1);
  OhmmsAttributeSet aAttrib;
  aAttrib.add(norb,"orbitals");
  aAttrib.add(degeneracy,"degeneracy");
  aAttrib.add(hrefname,"href");
  aAttrib.put(cur);
  if(norb ==0)
  {
    app_error() << "SplineSetBuilder::createSPOSet failed. Check the attribte orbitals." << std::endl;
    return 0;
  }
  app_log() << "    Degeneracy = " << degeneracy << std::endl;
  std::vector<int> npts(3);
  npts[0]=GridXYZ->nX;
  npts[1]=GridXYZ->nY;
  npts[2]=GridXYZ->nZ;
  std::vector<RealType> inData(npts[0]*npts[1]*npts[2]);
  SPOSetType* psi= new SPOSetType(norb);
  std::vector<int> occSet(norb);
  for(int i=0; i<norb; i++)
    occSet[i]=i;
  cur=cur->children;
  while(cur != NULL)
  {
    std::string cname((const char*)(cur->name));
    if(cname == "occupation")
    {
      std::string occ_mode("ground");
      const xmlChar* o=xmlGetProp(cur,(const xmlChar*)"mode");
      if(o!= NULL)
        occ_mode = (const char*)o;
      //Do nothing if mode == ground
      if(occ_mode == "excited")
      {
        std::vector<int> occ_in, occRemoved;
        putContent(occ_in,cur);
        for(int k=0; k<occ_in.size(); k++)
        {
          if(occ_in[k]<0)
            occRemoved.push_back(-occ_in[k]-1);
        }
        int kpopd=0;
        for(int k=0; k<occ_in.size(); k++)
        {
          if(occ_in[k]>0)
            occSet[occRemoved[kpopd++]]=occ_in[k]-1;
        }
      }
      hid_t h_file = H5Fopen(hrefname.c_str(),H5F_ACC_RDWR,H5P_DEFAULT);
      const xmlChar* h5path = xmlGetProp(cur,(const xmlChar*)"h5path");
      std::string hroot("/eigenstates_3/twist_0");
      if(h5path != NULL)
        hroot=(const char*)h5path;
      char wfname[128],wfshortname[16];
      for(int iorb=0; iorb<norb; iorb++)
      {
        sprintf(wfname,"%s/band_%d/eigenvector",hroot.c_str(),occSet[iorb]/degeneracy);
        sprintf(wfshortname,"b%d",occSet[iorb]/degeneracy);
        SPOType* neworb=0;
        std::map<std::string,SPOType*>::iterator it(NumericalOrbitals.find(wfshortname));
        if(it == NumericalOrbitals.end())
        {
          neworb=new SPOType(GridXYZ);
          HDFAttribIO<std::vector<RealType> > dummy(inData,npts);
          dummy.read(h_file,wfname);
          //neworb->reset(inData.begin(), inData.end(), targetPtcl.Lattice.BoxBConds[0]);
          neworb->reset(inData.begin(), inData.end(), targetPtcl.Lattice.SuperCellEnum);
          NumericalOrbitals[wfshortname]=neworb;
          app_log() << "   Reading spline function " << wfname << std::endl;
        }
        else
        {
          neworb = (*it).second;
          app_log() << "   Reusing spline function " << wfname << std::endl;
        }
        psi->add(neworb);
      }
      H5Fclose(h_file);
    }
    cur=cur->next;
  }
  SPOType* aorb=(*NumericalOrbitals.begin()).second;
  std::string fname("spline3d.vti");
  std::ofstream dfile(fname.c_str());
  dfile.setf(std::ios::scientific, std::ios::floatfield);
  dfile.setf(std::ios::left,std::ios::adjustfield);
  dfile.precision(10);
  dfile << "<?xml version=\"1.0\"?>" << std::endl;
  dfile << "<VTKFile type=\"ImageData\" version=\"0.1\">" << std::endl;
  dfile << "  <ImageData WholeExtent=\"0 " << npts[0]-2 << " 0 " << npts[1]-2 << " 0 " << npts[2]-2
        << "\" Origin=\"0 0 0\" Spacing=\"1 1 1\">"<< std::endl;
  dfile << "    <Piece Extent=\"0 " << npts[0]-2 << " 0 " << npts[1]-2 << " 0 " << npts[2]-2 << "\">" << std::endl;
  dfile << "       <PointData Scalars=\"wfs\">" << std::endl;
  dfile << "          <DataArray type=\"Float32\" Name=\"wfs\">" << std::endl;
  int ng=0;
  GradType grad;
  ValueType lap;
  for(int ix=0; ix<npts[0]-1; ix++)
  {
    double x(GridXYZ->gridX->operator()(ix));
    for(int iy=0; iy<npts[1]-1; iy++)
    {
      double y(GridXYZ->gridY->operator()(iy));
      for(int iz=0; iz<npts[2]-1; iz++, ng++)
      {
        PosType p(x,y,GridXYZ->gridZ->operator()(iz));
        //aorb.setgrid(p);
        //Timing with the std::ofstream is not correct.
        //Uncomment the line below and comment out the next two line.
        //double t=aorb.evaluate(p,grad,lap);
        dfile << std::setw(20) << aorb->evaluate(p,grad,lap);
        if(ng%5 == 4)
          dfile << std::endl;
      }
    }
  }
  dfile << "          </DataArray>" << std::endl;
  dfile << "       </PointData>" << std::endl;
  dfile << "    </Piece>" << std::endl;
  dfile << "  </ImageData>" << std::endl;
  dfile << "</VTKFile>" << std::endl;
  abort();
  return psi;
}
}
