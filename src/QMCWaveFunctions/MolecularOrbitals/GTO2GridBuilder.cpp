//////////////////////////////////////////////////////////////////
// (c) Copyright 2003- by Jeongnim Kim
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
#include "Numerics/Transform2GridFunctor.h"
#include "Numerics/OneDimCubicSpline.h"
#include "Numerics/HDFNumericAttrib.h"
#include "Numerics/GaussianBasisSet.h"
#include "QMCWaveFunctions/MolecularOrbitals/GridMolecularOrbitals.h"
#include "QMCWaveFunctions/MolecularOrbitals/GTO2GridBuilder.h"

namespace ohmmsqmc {

  bool
  GTO2GridBuilder::addRadialOrbital(xmlNodePtr cur, 
				    const QuantumNumberType& nlms) {

    int n=nlms[0];
    int l=nlms[1];

    //Using default <radfunc exponent="alpha" contraction="c"/>
    GaussianCombo<ValueType> gaussian(l,Normalized);
    gaussian.putBasisGroup(cur);

    //pointer to the grid
    GridType* agrid = m_orbitals->Grids[0];
    RadialOrbitalType *radorb = new OneDimCubicSpline<ValueType>(agrid);

    //spline the slater type orbital
    Transform2GridFunctor<GaussianCombo<RealType>,RadialOrbitalType> transform(gaussian, *radorb);
    transform.generate(agrid->rmin(), agrid->rmax(),agrid->size());

    //add the radial orbital to the list
    m_orbitals->Rnl.push_back(radorb);
    m_orbitals->RnlID.push_back(nlms);

    return true;
  }

  /** Default function to add a radial grid to the list of radial grids.
   * \param cur the current xmlNode to be processed
   * \return true if succeeds
   */
  bool 
  GTO2GridBuilder::addGrid(xmlNodePtr cur) {

    if(!m_orbitals) {
      ERRORMSG("m_orbitals, SphericalOrbitals<ROT,GT>*, is not initialized")
      return false;
    }

    XMLReport("Converting analytic orbitals to radial grid functions. Modify to use zero-based grid.")
    RealType ri = 1e-6;
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
