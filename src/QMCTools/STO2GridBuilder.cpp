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
#include "QMCWaveFunctions/MolecularOrbitals/GridMolecularOrbitals.h"
#include "QMCWaveFunctions/MolecularOrbitals/STO2GridBuilder.h"
using namespace std;
namespace qmcplusplus {

  bool
  STO2GridBuilder::putCommon(xmlNodePtr cur) {
    return true;
  }

  /** Add a new Slater Type Orbital with quantum numbers \f$(n,l,m,s)\f$ 
   * \param cur  the current xmlNode to be processed
   * \param nlms a vector containing the quantum numbers \f$(n,l,m,s)\f$
   * \return true is succeeds
   *
   This function puts the STO on a logarithmic grid and calculates the boundary 
   conditions for the 1D Cubic Spline.  The derivates at the endpoint 
   are assumed to be all zero.  Note: for the radial orbital we use
   \f[ f(r) = \frac{R(r)}{r^l}, \f] where \f$ R(r) \f$ is the usual
   radial orbital and \f$ l \f$ is the angular momentum.
  */
  bool
  STO2GridBuilder::addRadialOrbital(xmlNodePtr cur, 
				    const QuantumNumberType& nlms) {
    if(!m_orbitals) {
      ERRORMSG("m_orbitals, SphericalOrbitals<ROT,GT>*, is not initialized")
      return false;
    }

    RadialOrbitalType *radorb =  NULL;
    int n=nlms[0];
    int l=nlms[1];
    RealType zeta = 1.0;
    xmlNodePtr s = cur->xmlChildrenNode;
    while(s != NULL) {
      string cname((const char*)(s->name));
      if(cname == "parameter" || cname =="Var") {
	putContent(zeta,s);
      }
      s=s->next;
    }
    XMLReport("Zeta = " << zeta)
    STONorm<RealType> anorm(n);
    
    GenericSTO<RealType> sto(n-l-1,zeta,anorm(n-1,zeta));
    XMLReport("Calculating 1D-Cubic spline.")

      //pointer to the grid
    GridType* agrid = m_orbitals->Grids[0];

    radorb = new OneDimCubicSpline<RealType>(agrid);
    //spline the slater type orbital
    Transform2GridFunctor<GenericSTO<RealType>,RadialOrbitalType> transform(sto, *radorb);
    transform.generate(agrid->rmin(), agrid->rmax(),agrid->size());
    //add the radial orbital to the list
    m_orbitals->Rnl.push_back(radorb);
    m_orbitals->RnlID.push_back(nlms);

    return true;
  }
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
