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
#include "Numerics/SlaterTypeOrbital.h"
#include "Numerics/Transform2GridFunctor.h"
#include "Numerics/OneDimCubicSpline.h"
#include "Numerics/HDFNumericAttrib.h"
#include "QMCWaveFunctions/MolecularOrbitals/GridMolecularOrbitals.h"
#include "QMCWaveFunctions/MolecularOrbitals/RGFBuilderBase.h"
using namespace std;
namespace ohmmsqmc {

  /** the destructor
  */
  RGFBuilderBase::~RGFBuilderBase() {
    //clean up
    H5Gclose(m_group_id);   
    H5Fclose(m_file_id);   
  }

  void RGFBuilderBase::print(const string& aroot, int omode) {

    enum {DEBUG_OUTPUT=0, ASCII_OUTPUT, HDF_OUTPUT};

    if(omode == DEBUG_OUTPUT) {
       GridType* agrid = m_orbitals->Grids[0];
      for(int i=0; i<m_orbitals->Rnl.size(); i++) {
        RadialOrbitalType* radorb = m_orbitals->Rnl[i];
        char fname[128];
        sprintf(fname,"%s.%d.dat",aroot.c_str(), i);
        ofstream dfile(fname);
        dfile.setf(ios::scientific, ios::floatfield);
        dfile.setf(ios::left,ios::adjustfield);
        dfile.precision(12);
        for(int ig=0; ig<agrid->size()-1; ig++) {
          RealType dr = ((*agrid)(ig+1)-(*agrid)(ig))/5.0;
          RealType _r =(*agrid)(ig)+dr*0.00131,y,dy,d2y;
          while(_r<(*agrid)(ig+1)) {
            radorb->setgrid(_r);
            y = radorb->evaluate(_r,1.0/_r,dy,d2y);
            dfile << setw(24) << _r << setw(24) << y 
                  << setw(24) << dy << setw(24) << d2y
                  << endl;
            _r+=dr;
          }
        }
      }
    } else if(omode == ASCII_OUTPUT) {
      GridType* agrid = m_orbitals->Grids[0];
      char fname[128];
      sprintf(fname,"%s.combo.dat",aroot.c_str());
      ofstream dfile(fname);
      dfile.setf(ios::scientific, ios::floatfield);
      dfile.setf(ios::left,ios::adjustfield);
      int norb = m_orbitals->Rnl.size();
      dfile << "#column 0 : Radial grid" << endl;
      for(int i=0; i<norb; i++) {
        dfile << "#column " << i+1 << " :  " << m_orbitals->RnlID[i] << endl;
      }
      dfile.precision(15);
      for(int ig=0; ig<agrid->size()-1; ig++) {
        RealType _r = (*agrid)(ig);
        RealType _rinv = 1.0/_r;
        agrid->index(_r);
        dfile << setw(30) << _r;
        for(int i=0; i<norb; i++) {
          dfile << setw(30) << m_orbitals->Rnl[i]->evaluate(_r,_rinv);
        }
        dfile << endl;
      }
    }
  }

  /** Default function to add a radial grid to the list of radial grids.
   * \param cur the current xmlNode to be processed
   * \return true if succeeds
   */
  bool 
  RGFBuilderBase::addGrid(xmlNodePtr cur) {
    if(!m_orbitals) {
      ERRORMSG("m_orbitals, SphericalOrbitals<ROT,GT>*, is not initialized")
      return false;
    }

    XMLReport("Converting analytic orbitals to radial grid functions")
    RealType ri = 1e-5;
    RealType rf = 100.0;
    IndexType npts = 1001;
    const xmlChar* ri_ptr = xmlGetProp(cur,(const xmlChar *)"ri");
    const xmlChar* rf_ptr = xmlGetProp(cur,(const xmlChar *)"rf");
    const xmlChar* n_ptr = xmlGetProp(cur,(const xmlChar *)"npts");

    if(ri_ptr) ri = atof((const char*)ri_ptr);
    if(rf_ptr) rf = atof((const char*)rf_ptr);
    if(n_ptr) npts = atoi((const char*)n_ptr);
    LOGMSG("Using log grid with default values: ri = " << ri << " rf = " << rf << " npts = " << npts)

    GridType *agrid = new LogGrid<RealType>;
    agrid->set(ri,rf,npts);

    m_orbitals->Grids.push_back(agrid);
    return true;
  }
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
