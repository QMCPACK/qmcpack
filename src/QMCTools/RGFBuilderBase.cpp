//////////////////////////////////////////////////////////////////
// (c) Copyright 2003 by Jeongnim Kim and Jordan Vincent
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
#include "Utilities/OhmmsInfo.h"
#include "Numerics/LibxmlNumericIO.h"
#include "Numerics/OneDimCubicSpline.h"
#if defined(HAVE_LIBHDF5)
#include "Numerics/HDFNumericAttrib.h"
#endif
#include "QMCTools/GridMolecularOrbitals.h"
#include "QMCTools/RGFBuilderBase.h"
#include "QMCFactory/OneDimGridFactory.h"
using namespace std;
namespace qmcplusplus
{

/** the destructor
*/
RGFBuilderBase::~RGFBuilderBase()
{
#if defined(HAVE_LIBHDF5)
  if(m_group_id>-1)
    H5Gclose(m_group_id);
  if(m_file_id>-1)
    H5Fclose(m_file_id);
#endif
}

/** printout for debug
 * @param aroot filename
 * @param omode mode, 0=check cubic spline quality, 1=ascii and 2=hdf5
 */
void RGFBuilderBase::print(const string& aroot, int omode)
{
  enum {DEBUG_OUTPUT=0, ASCII_OUTPUT, HDF_OUTPUT};
  if(omode == DEBUG_OUTPUT)
  {
    GridType& agrid(m_orbitals->Rnl[0]->grid());
    for(int i=0; i<m_orbitals->Rnl.size(); i++)
    {
      RadialOrbitalType* radorb = m_orbitals->Rnl[i];
      char fname[128];
      sprintf(fname,"%s.%d.dat",aroot.c_str(), i);
      ofstream dfile(fname);
      dfile.setf(ios::scientific, ios::floatfield);
      dfile.setf(ios::left,ios::adjustfield);
      dfile.precision(12);
      for(int ig=0; ig<agrid.size()-1; ig++)
      {
        RealType dr = (agrid(ig+1)-agrid(ig))/5.0;
        RealType _r =agrid(ig)+dr*0.00131,y,dy,d2y;
        while(_r<agrid(ig+1))
        {
          //setgrid is removed.  The GridManager functor will re-evaluate the grid parameters.
          //radorb->setgrid(_r);
          y = radorb->evaluate(_r,1.0/_r,dy,d2y);
          dfile << setw(24) << _r << setw(24) << y
                << setw(24) << dy << setw(24) << d2y
                << endl;
          _r+=dr;
        }
      }
    }
  }
  else
    if(omode == ASCII_OUTPUT)
    {
      GridType& agrid(m_orbitals->Rnl[0]->grid());
      char fname[128];
      sprintf(fname,"%s.combo.dat",aroot.c_str());
      ofstream dfile(fname);
      dfile.setf(ios::scientific, ios::floatfield);
      dfile.setf(ios::left,ios::adjustfield);
      int norb = m_orbitals->Rnl.size();
      dfile << "# column 0 : Radial grid" << endl;
      for(int i=0; i<norb; i++)
      {
        dfile << "# column " << i+1 << " :  " << m_orbitals->RnlID[i] << endl;
      }
      dfile.precision(15);
      for(int ig=0; ig<agrid.size()-1; ig++)
      {
        RealType _r = agrid(ig);
        dfile << setw(30) << _r;
        for(int i=0; i<norb; i++)
        {
          dfile << setw(30) << m_orbitals->Rnl[i]->operator()(ig);
        }
        dfile << endl;
      }
    }
}

/** Add a radial grid to the list of radial grids.
 * \param cur current xmlNode to be processed
 * \return true if succeeds
 *
 * The default grid is LogGrid and  the number of grids is 1001.
 * Valid radial grids instantiated by the attributes
 * - LogGrid
 *  -- \<grid type="log" ri="first-grid-value" rf="final-grid-value" npts="grid-points"/\>
 *  -- \f$r(i) = ri\times(rf/ri)^{i/(npts-1)}\f$ for i=[0,npts).
 *  -- Default values : ri=1e-5, rf=100
 * - LogGridZero
 *  -- \<grid type="log" scale="scale-value" step="step-value" first="1" npts="grid-points"/\>
 *  -- \f$r(i) = scale\times(\exp{step*i}-1.0)\f$ for i=[0,npts).
 *  -- Follows Pseudopotential xml schema (fsatom-pp)
 *  -- All the attributes should be provided.first attribute is ignored.
 * - LinearGrid
 *  -- \<grid type="linear" ri="first-grid-value" rf="final-grid-value" npts="grid-points"/\>
 */
bool
RGFBuilderBase::addGrid(xmlNodePtr cur)
{
  if(!m_orbitals)
  {
    ERRORMSG("m_orbitals, SphericalOrbitals<ROT,GT>*, is not initialized")
    return false;
  }
  XMLReport("Converting analytic orbitals to radial grid functions")
  //GridType *agrid=0;
  //RealType ri = 1e-5;
  //RealType rf = 100.0;
  //RealType ascale = -1.0e0;
  //RealType astep = 1.25e-2;
  //IndexType npts = 1001;
  //string gridType("log");
  //OhmmsAttributeSet radAttrib;
  //radAttrib.add(gridType,"type");
  //radAttrib.add(npts,"npts");
  //radAttrib.add(ri,"ri"); radAttrib.add(rf,"rf");
  //radAttrib.add(ascale,"ascale"); radAttrib.add(astep,"astep");
  //if(cur) radAttrib.put(cur);
  //if(gridType == "log") {
  //  if(ascale>0.0) {
  //    LOGMSG("Using log grid with default values: scale = " << ascale << " step = " << astep << " npts = " << npts)
  //    agrid = new LogGridZero<RealType>;
  //    agrid->set(astep,ascale,npts);
  //  } else {
  //    LOGMSG("Using log grid with default values: ri = " << ri << " rf = " << rf << " npts = " << npts)
  //    agrid = new LogGrid<RealType>;
  //    agrid->set(ri,rf,npts);
  //  }
  //} else if(gridType == "linear") {
  //  LOGMSG("Using linear grid with default values: ri = " << ri << " rf = " << rf << " npts = " << npts)
  //  agrid = new LinearGrid<RealType>;
  //  agrid->set(ri,rf,npts);
  //}
  GridType *agrid = OneDimGridFactory::createGrid(cur);
  m_orbitals->Grids.push_back(agrid);
  return true;
}
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/
