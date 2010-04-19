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
#if QMC_BUILD_LEVEL>2 && OHMMS_DIM==3
#include "QMCWaveFunctions/Fermion/BackflowTransformation.h"
#endif
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
    BasisSetFactory* myBasisSetFactory;
    SlaterDeterminant_t* slaterdet_0;
    
#if QMC_BUILD_LEVEL>2 && OHMMS_DIM==3
    bool UseBackflow;
    BackflowTransformation *BFTrans;
#endif

    /** process a determinant element
     * @param cur xml node
     * @param firstIndex index of the determinant
     * @return firstIndex+number of orbitals
     */
    bool putDeterminant(xmlNodePtr cur, int firstIndex);
  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
