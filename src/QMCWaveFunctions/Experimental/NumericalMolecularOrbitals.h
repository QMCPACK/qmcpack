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
#ifndef QMCPLUSPLUS_NUMERICAL_MOLECULARORBITALS_H
#define QMCPLUSPLUS_NUMERICAL_MOLECULARORBITALS_H
#include "QMCWaveFunctions/OrbitalBuilderBase.h"

namespace qmcplusplus {

  class GridMolecularOrbitals;

  class NumericalMolecularOrbitals: public OrbitalBuilderBase {

    GridMolecularOrbitals *Original;

    public:

    /** constructor
     * \param wfs reference to the wavefunction
     * \param ions reference to the ions
     * \param els reference to the electrons
     */
    NumericalMolecularOrbitals(ParticleSet& els, TrialWaveFunction& wfs, ParticleSet& ions);

    /** initialize the Antisymmetric wave function for electrons
     *@param cur the current xml node
     *
     */
    bool put(xmlNodePtr cur);

  private:

  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
