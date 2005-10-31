//////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Jeongnim Kim
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
#ifndef QMCPLUSPLUS_NUMERICALORBITALSETBUILDER_H
#define QMCPLUSPLUS_NUMERICALORBITALSETBUILDER_H

#include "Utilities/OhmmsInfo.h"
#include "QMCWaveFunctions/OrbitalBuilderBase.h"
#include "Numerics/TriCubicSplineT.h"
#include "QMCWaveFunctions/MixedSPOSet.h"
#include "QMCWaveFunctions/LCOrbitals.h"
#include "QMCWaveFunctions/MolecularOrbitals/GridMolecularOrbitals.h"
#include "QMCWaveFunctions/DiracDeterminant.h"
#include "QMCWaveFunctions/SlaterDeterminant.h"

namespace qmcplusplus {

  /**@ingroup WFSBuilder
   * A builder class for a set of SlaterDeterminant with mixed orbitals
   *
   * The most general Orbital Builder for typical electronic
   * structure calculations with ions.
   */
  class NumericalOrbitalSetBuilder: public OrbitalBuilderBase {

  public:

    //@typedef M(olecular) O(rbital) Basis Set Type
    typedef GridMolecularOrbitals::BasisSetType MOBasisSetType;
    //@typedef L(ocalized) O(rbital) Type
    typedef LCOrbitals<MOBasisSetType> LOType;
    //@typedef S(ingle) P(article) O(rbital) Set Type
    typedef MixedSPOSet<LOType> SPOSetType;
    //@typedef N(umerical) O(rbital) Type
    typedef TriCubicSplineT<ValueType,RealType> NumericalOrbitalType;
    //@typedef Dirac Determinant Type using MixedSPOSet
    typedef DiracDeterminant<SPOSetType>  Det_t;
    //@typedef Slater Determinant Type using MixedSPOSet
    typedef SlaterDeterminant<SPOSetType> SlaterDeterminant_t;

    /** constructor
     * @param p target ParticleSet
     * @param psi TrialWaveFunction 
     * @param psets a set of ParticleSet objects
     */
    NumericalOrbitalSetBuilder(ParticleSet& p, 
        TrialWaveFunction& psi, PtclPoolType& psets);

    /** initialize the Antisymmetric wave function for electrons
     *@param cur the current xml node
     */
    bool put(xmlNodePtr cur);

  private:

    ///handle localized basis set
    GridMolecularOrbitals*  LocalizedBasisBuilder;
    ///reference to a ParticleSetPool
    PtclPoolType& ptclPool;
    ///a global grid for all the Numerical Orbitals that are created by this builder
    XYZCubicGrid<RealType>* GridXYZ;
    ///Molecular Orbital Basis Set
    MOBasisSetType* MOBasisSet;
    /** set of numerical orbitals 
     *
     * NGOrbitals[name] points to a unique numerical orbital.
     */
    map<string,NumericalOrbitalType* > NumericalOrbitals;
    /** set of localized orbitals 
     */
    map<string,LOType*>                LocalizedOrbitals;

    /** set of Single-particle orbital set for
     */
    map<string,SPOSetType*>       SPOSet;

    /** set of Diract Determinants */
    map<string,Det_t*>            DetSet;

    /** set of SlaterDeterminant for multi determinant cases
     */
    vector<SlaterDeterminant_t*>  SlaterDetSet;

    SlaterDeterminant_t* addSlaterDeterminant(xmlNodePtr cur);
    SPOSetType* createSPOSet(xmlNodePtr cur, const string& ref, int norb);
  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
