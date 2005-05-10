//////////////////////////////////////////////////////////////////
// (c) Copyright 2003  by Jeongnim Kim
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
// -*- C++ -*-
#ifndef OHMMS_QMC_NUMERICALORBITALSETBUILDER_H
#define OHMMS_QMC_NUMERICALORBITALSETBUILDER_H

#include "QMCWaveFunctions/OrbitalBuilderBase.h"
#include "Numerics/TriCubicSplineT.h"

namespace ohmmsqmc {

 /** builder for a set of SlaterDeterminant with numerical orbitals
  */
  class NumericalOrbitalSetBuilder: public OrbitalBuilderBase {

  public:

    ///typedef for Single-Particle-Orbita s
    typedef TriCubicSplineT<ValueType,RealType> SPOType;
    ///typedef for the set of Single-Particle-Orbitals
    typedef SingleParticleOrbitalSet<SPOType>   SPOSetType;

    /** constructor
     * \param wfs reference to the wavefunction
     * \param ions reference to the ions
     * \param els reference to the electrons
     */
    NumericalOrbitalSetBuilder(TrialWaveFunction& wfs);

    /** initialize the Antisymmetric wave function for electrons
     *@param cur the current xml node
     */
    bool put(xmlNodePtr cur);

  private:

    ///a global grid for all the Numerical Orbitals that are created by this builder
    XYZCubicGrid<RealType> *GridXYZ;

    /** a set of numerical orbitals 
     *
     * SPOSet[name] points to a unique numerical orbital.
     */
    map<string,TriCubicSplineT<ValueType>* > SPOSet;

    bool addSlaterDeterminantSet(cur);
  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
