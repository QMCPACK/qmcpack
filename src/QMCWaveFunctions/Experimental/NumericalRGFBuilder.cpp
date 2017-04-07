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
    
    
#include "Utilities/OhmmsInfo.h"
#include "Numerics/LibxmlNumericIO.h"
#include "Numerics/SlaterTypeOrbital.h"
#include "Numerics/Transform2GridFunctor.h"
#include "Numerics/OneDimCubicSpline.h"
#include "QMCWaveFunctions/MolecularOrbitals/GridMolecularOrbitals.h"
#include "QMCWaveFunctions/MolecularOrbitals/NumericalRGFBuilder.h"
#if defined(HAVE_LIBHDF5)
#include "Numerics/HDFNumericAttrib.h"
#endif
#include "OhmmsData/AttributeSet.h"

namespace qmcplusplus
{

NumericalRGFBuilder::NumericalRGFBuilder(xmlNodePtr cur)
{
  if(cur != NULL)
    putCommon(cur);
}

bool NumericalRGFBuilder::putCommon(xmlNodePtr cur)
{
#if defined(HAVE_LIBHDF5)
  std::string afilename("invalid");
  OhmmsAttributeSet aAttrib;
  aAttrib.add(afilename,"href");
  aAttrib.add(afilename,"file");
  aAttrib.put(cur);
  if(afilename == "invalid")
  {
    ERRORMSG("NumericalRGFBuilder::NumericalRGFBuilder: Missing a file name")
  }
  else
  {
    //open the HDF5 file
    m_file_id = H5Fopen(afilename.c_str(),H5F_ACC_RDWR,H5P_DEFAULT);
    XMLReport("Opening file: " << afilename)
  }
#endif
  return true;
}

/**  Add a new numerical radial orbital with quantum numbers \f$(n,l,m,s)\f$ to the list of radial orbitals.
* \param cur the current xmlNode to be processed
* \param nlms a vector containing the quantum numbers \f$(n,l,m,s)\f$
* \return true is succeeds
*
This function finds the location of the corresponding data for
the radial orbital in the HDF5 file, reads in the data into a new
RadialOrbitalType and calculates the boundary conditions for the
1D Cubic Spline.  The derivates at the endpoint are assumed to be
all zero.  Note: for the radial orbital we use
\f[ f(r) = \frac{R(r)}{r^l}, \f] where \f$ R(r) \f$ is the usual
radial orbital and \f$ l \f$ is the angular momentum.
*/
bool
NumericalRGFBuilder::addRadialOrbital(xmlNodePtr cur,
                                      const QuantumNumberType& nlms)
{
  if(!m_orbitals)
  {
    ERRORMSG("m_orbitals, SphericalOrbitals<ROT,GT>*, is not initialized")
    return false;
  }
  int n=nlms[0];
  int l=nlms[1];
  int m=nlms[2];
  int s=nlms[3];
  //first valid index for the radial grid
  int imin = 0;
  xmlAttrPtr att = cur->properties;
  while(att != NULL)
  {
    std::string aname((const char*)(att->name));
    if(aname == "imin")
    {
      imin = atoi((const char*)(att->children->content));
      XMLReport("First valid index for radial grid = " << imin)
    }
    att = att->next;
  }
#if defined(HAVE_LIBHDF5)
  //open the group containing the proper orbital
  char grpname[128];
  sprintf(grpname,"orbital%04d",Counter++);
  hid_t group_id_orb = H5Gopen(m_group_id,grpname);
  /*find the correct power such that:
    R(r)/r^l = rad_orb(r)/r^{power}*/
  Vector<int> Power;
  HDFAttribIO<Vector<int> > PowerHDFIn(Power);
  PowerHDFIn.read(group_id_orb,"power");
  //check that the quantum numbers are consistent
  Vector<int> quantum_numbers;
  HDFAttribIO<Vector<int> > QuantumNoHDFIn(quantum_numbers);
  QuantumNoHDFIn.read(group_id_orb,"quantum_numbers");
  if((n != quantum_numbers[0]) || (l != quantum_numbers[1]))
  {
    ERRORMSG("Inconsistency between input file and HDF5 file!")
    return false;
  }
  XMLReport("Reading from HDF5 file for Zeta = " << quantum_numbers[2])
  //read the radial orbital data
  Vector<double> rad_orb;
  HDFAttribIO<Vector<double> > RadOrbIn(rad_orb);
  RadOrbIn.read(group_id_orb,"radial_orbital");
  H5Gclose(group_id_orb);
  //pointer to the grid
  GridType* agrid = m_orbitals->Grids[0];
  int rinv_p = Power[0];
  if(rinv_p != 0)
    for(int ig=0; ig<rad_orb.size(); ig++)
      rad_orb[ig] *= std::pow(agrid->r(ig),-rinv_p);
  XMLReport("Multiplying orbital by r^" << -rinv_p)
  //last valid index for radial grid
  int imax = rad_orb.size()-1;
  RadialOrbitalType *radorb = new OneDimCubicSpline<RealType>(agrid,rad_orb);
  //calculate boundary condition, assume derivates at endpoint are 0.0
  RealType yprime_i = rad_orb[imin+1]-rad_orb[imin];
  if(std::abs(yprime_i)<1e-10)
    yprime_i = 0.0;
  yprime_i /= (agrid->r(imin+1)-agrid->r(imin));
  //set up 1D-Cubic Spline
  radorb->spline(imin,yprime_i,imax,0.0);
  int lastRnl=m_orbitals->Rnl.size();
  //add radial orbital to list
  m_orbitals->Rnl.push_back(radorb);
  m_orbitals->RnlID.push_back(nlms);
  XMLReport("Calculating 1D-Cubic spline, cusp condtion = " << yprime_i)
  //#ifdef PRINT_DEBUG
  std::string fname(grpname);
  fname.append(".dat");
  std::ofstream dfile(fname.c_str());
  for(int ig=imin; ig<radorb->size(); ig++)
  {
    RealType dr = (radorb->r(ig+1)- radorb->r(ig))/5.0;
    RealType _r = radorb->r(ig),y,dy,d2y;
    while(_r<radorb->r(ig+1))
    {
      //Do not need this anymore
      //radorb->setgrid(_r);
      y = radorb->evaluate(_r,1.0/_r,dy,d2y);
      dfile << std::setw(15) << _r << std::setw(20) << setprecision(12) << y
            << std::setw(20) << dy << std::setw(20) << d2y
            << std::endl;
      _r+=dr;
    }
  }
  //#endif
  //
  if(lastRnl)
  {
    XMLReport("\tSetting the GridManager of " << lastRnl << " radial orbital to false")
    radorb->setGridManager(false);
  }
  return true;
#else
  ERRORMSG("HDF5 is disabled.")
  OHMMS::Controller->abort();
#endif
}

/** Add a radial grid to the list of radial grids from a HDF5 file
  * \param cur the current xmlNode to be processed
  * \return true if succeeds
  *
  This function opens the HDF5 file named in the input file,
  finds the location of the radial grid data and creates the
  grid.  The possible grid types are:
  Linear \f[ r_i = r_0 +
  i\left( \frac{r_f - r_0}{N-1} \right) \f]
  Logarithmic \f[ r_i = r_0
  \left( \frac{r_f}{r_0} \right) ^{\frac{i}{N-1}} \f]
  or LogZero \f[ r_i = B
  \left[ \exp(Ai)-1 \right] \f]
  @todo use OneDimGridFactory
*/
bool
NumericalRGFBuilder::addGrid(xmlNodePtr cur)
{
  if(!m_orbitals)
  {
    ERRORMSG("m_orbitals, SphericalOrbitals<ROT,GT>*, is not initialized")
    return false;
  }
#if defined(HAVE_LIBHDF5)
  //create a grid and initialize it
  GridType* agrid = 0;
  if(m_file_id>=0)
  {
    //const xmlChar* fptr = xmlGetProp(cur, (const xmlChar *)"file");
    //if(fptr) {
    //open the HDF5 file
    //string afilename((const char*)fptr);
    //m_file_id = H5Fopen(afilename.c_str(),H5F_ACC_RDWR,H5P_DEFAULT);
    //XMLReport("Opening file: " << afilename)
    m_group_id = H5Gopen(m_file_id,"radial_basis_states");
    hid_t group_id_grid = H5Gopen(m_group_id,"grid");
    Vector<int> grid_type;
    Vector<double> grid_params;
    HDFAttribIO<Vector<int> > GridTypeIn(grid_type);
    GridTypeIn.read(group_id_grid,"grid_type");
    HDFAttribIO<Vector<double> > GridParamsIn(grid_params);
    GridParamsIn.read(group_id_grid,"params");
    /* documentation */
    enum {LINEAR_GRID=1, LOG_GRID, LOG_GRID_ZERO};
    H5Gclose(group_id_grid);
    if(grid_type[0] == LINEAR_GRID)
    {
      agrid = new LinearGrid<RealType>;
      LOGMSG("Using a linear grid.")
      RealType ri = grid_params(0);
      RealType rf = grid_params(1);
      IndexType npts = static_cast<int>(grid_params(2));
      agrid->set(ri,rf,npts);
    }
    else
      if(grid_type[0] == LOG_GRID)
      {
        agrid = new LogGrid<RealType>;
        LOGMSG("Using a log grid.")
        RealType ri = grid_params(0);
        RealType rf = grid_params(1);
        IndexType npts = static_cast<int>(grid_params(2));
        agrid->set(ri,rf,npts);
      }
      else
        if(grid_type[0] == LOG_GRID_ZERO)
        {
          agrid = new LogGridZero<RealType>;
          LOGMSG("Using a log zero grid.")
          RealType a = grid_params(0);
          RealType b = grid_params(1);
          IndexType npts = static_cast<int>(grid_params(2));
          agrid->set(a,b,npts);
        }
    LOGMSG("Grid Values: ri = " << agrid->rmin() << " rf = " << agrid->rmax() << " npts = " << agrid->size())
    m_orbitals->Grids.push_back(agrid);
  }
  //exit(-1);
  return true;
#else
  ERRORMSG("HDF5 is disabled.")
  OHMMS::Controller->abort();
#endif
}

}
