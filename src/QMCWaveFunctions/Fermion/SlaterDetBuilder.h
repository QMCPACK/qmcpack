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
#ifndef QMCPLUSPLUS_LCORBITALSETBUILDER_H
#define QMCPLUSPLUS_LCORBITALSETBUILDER_H

#include "QMCWaveFunctions/OrbitalBuilderBase.h"
#include "QMCWaveFunctions/BasisSetFactory.h"
#include "QMCWaveFunctions/Fermion/SlaterDet.h"

namespace qmcplusplus {

 /** derived class from OrbitalBuilderBase
  *
  * Builder SlaterDeterminant with LCOrbitalSet
  */
  class SlaterDetBuilder: public OrbitalBuilderBase {

  public:

    typedef SlaterDet SlaterDeterminant_t;
    typedef DiracDeterminantBase Det_t;
    /** constructor
     * \param els reference to the electrons
     * \param psi reference to the wavefunction
     * \param ions reference to the ions
     */
    SlaterDetBuilder(ParticleSet& els, TrialWaveFunction& psi, PtclPoolType& psets);

    ~SlaterDetBuilder();

    /** initialize the Antisymmetric wave function for electrons
     *@param cur the current xml node
     *
     */
    bool put(xmlNodePtr cur);

  private:

    ///reference to a PtclPoolType
    PtclPoolType& ptclPool;

    ////map<string,BasisSetBase*> BasisSet;
    //BasisSetBase* myBasisSet;

    map<string,SPOSetBasePtr> SPOSet;

    map<string,Det_t*> DetSet;

    vector<SlaterDeterminant_t*> SlaterDetSet;

    vector<RealType> sdet_coeff;

    /** process a determinant element
     * @param cur xml node
     * @param firstIndex index of the determinant
     * @return firstIndex+number of orbitals
     */
    int putDeterminant(xmlNodePtr cur, int firstIndex);

    /** build a Slater Determinant
     */
    void buildSlaterDetermiant();

    /** build a Multi-Slater Determinant
     */
    void buildMultiSlaterDetermiant();

    BasisSetFactory* myBasisSetFactory;
  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
