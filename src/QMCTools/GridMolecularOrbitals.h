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
#ifndef QMCPLUSPLUS_CONVERT2_RADIALGRID_H
#define QMCPLUSPLUS_CONVERT2_RADIALGRID_H

#include "QMCWaveFunctions/OrbitalBuilderBase.h"
#include "QMCTools/MolecularOrbitalBasis.h"
#include "QMCTools/RGFBuilderBase.h"

namespace qmcplusplus {

 /** derived class from OrbitalBuilderBase
  *
  * Create a basis set of molecular orbital types as defined in MolecularOrbitalBasis
  * with radial wave functions on the radial grids.
  */
  class GridMolecularOrbitals: public OrbitalBuilderBase {

  public:

    //@typedef radial grid type
    typedef OneDimGridBase<RealType>                       GridType;
    //@typedef \f$R_{nl}\f$ radial functor type defined on a grid
    typedef OneDimGridFunctor<RealType>                    RadialOrbitalType;
    //@typedef centered orbital type which represents a set of \f$R_{nl}Y_{lm}\f$  for a center
    typedef SphericalOrbitalSet<RadialOrbitalType,GridType> CenteredOrbitalType;
    //@typedef molecuar orbital basis composed of multiple CenteredOrbitalType s
    typedef MolecularOrbitalBasis<CenteredOrbitalType>      BasisSetType;

    /** constructor
     * \param els reference to the electrons
     * \param psi reference to the wavefunction
     * \param ions reference to the ions
     */
    GridMolecularOrbitals(ParticleSet& els, TrialWaveFunction& psi, ParticleSet& ions);

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

    enum {DONOT_EXPAND=0, GAUSSIAN_EXPAND=1, NATURAL_EXPAND};

    ///reference to the ionic system with which the basis set are associated
    ParticleSet& IonSys;

    ///pointer to the BasisSet built by GridMolecularOrbitals
    BasisSetType*      BasisSet;

    ///distance table pointer
    DistanceTableData* d_table;

    ///Current RadiaaGridFunctorBuilder
    RGFBuilderBase* rbuilder; 
    ///map for the radial orbitals
    map<string,int>    RnlID;
    ///map for the centers
    map<string,int>    CenterID;

    ///map for (n,l,m,s) to its quantum number index
    map<string,int> nlms_id;

    ///append Ylm channels
    int expandYlm(const string& rnl, const QuantumNumberType& nlms, int num, 
                  CenteredOrbitalType* aos, xmlNodePtr cur1, 
                  int expandlm=DONOT_EXPAND);

  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
