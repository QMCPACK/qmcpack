//////////////////////////////////////////////////////////////////
// (c) Copyright 2005-  by Jeongnim Kim
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
#ifndef QMCPLUSPLUS_ANY2RADIALGRIDFUNCTOR_H
#define QMCPLUSPLUS_ANY2RADIALGRIDFUNCTOR_H

#include "QMCWaveFunctions/MolecularOrbitals/RGFBuilderBase.h"

namespace qmcplusplus {

  /**Class to convert SlaterTypeOrbital to a radial orbital on a log grid.
   *
   * For a center,
   *   - only one grid is used
   *   - any number of radial orbitals 
   */
  struct Any2GridBuilder: public RGFBuilderBase {

    ///constructor
    Any2GridBuilder(xmlNodePtr cur=NULL);

    ///implement the virtual function
    bool addRadialOrbital(xmlNodePtr cur, const QuantumNumberType& nlms);
    bool putCommon(xmlNodePtr cur);

    bool Normalized;
    RealType m_rcut;
    QuantumNumberType m_nlms;

    void addGaussian(xmlNodePtr cur);
    void addSlater(xmlNodePtr cur);
    void addPade(xmlNodePtr cur);
  };

}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
