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
#ifndef OHMMS_QMC_RADIALGRID_MOLECULARORBITALS_H
#define OHMMS_QMC_RADIALGRID_MOLECULARORBITALS_H

#include <vector>
#include "Numerics/OneDimGridFunctor.h"
#include "QMCWaveFunctions/SphericalOrbitalSet.h"
#include "QMCWaveFunctions/MolecularOrbitals/MolecularOrbitalBasis.h"
#include "QMCWaveFunctions/LCOrbitals.h"
#include "QMCWaveFunctions/OrbitalBuilderBase.h"

namespace ohmmsqmc {

 /** derived class from OrbitalBuilderBase
  *
  * Create a basis set of molecular orbital types as defined in MolecularOrbitalBasis
  * with radial wave functions on the radial grids.
  */
  class GridMolecularOrbitals: public OrbitalBuilderBase {

  public:

    typedef OneDimGridBase<ValueType>                       GridType;
    typedef OneDimGridFunctor<ValueType>                    RadialOrbitalType;
    typedef SphericalOrbitalSet<RadialOrbitalType,GridType> CenteredOrbitalType;
    typedef MolecularOrbitalBasis<CenteredOrbitalType>      BasisSetType;


    /** constructor
     * \param wfs reference to the wavefunction
     * \param ions reference to the ions
     * \param els reference to the electrons
     */
    GridMolecularOrbitals(TrialWaveFunction& wfs, 
			  ParticleSet& ions, 
			  ParticleSet& els);

    /** initialize the Antisymmetric wave function for electrons
     *@param cur the current xml node
     *
     */
    bool put(xmlNodePtr cur);

    /** process basis element to build the basis set
     *@param cur the current xml node
     *@return a pointer to the BasisSet 
     *
     *This member function is necessary in order to use SDSetBuilderWithBasisSet
     *to initialize fermion wavefunctions with slater  determinants.
     */
    BasisSetType* addBasisSet(xmlNodePtr cur);

  private:

    ///pointer to the BasisSet built by GridMolecularOrbitals
    BasisSetType*      BasisSet;

    ///distance table pointer
    DistanceTableData* d_table;
    ///map for the radial orbitals
    map<string,int>    RnlID;
    ///map for the centers
    map<string,int>    CenterID;

  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
