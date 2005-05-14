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
#ifndef OHMMS_QMC_GENERIC_MOLECULARORBITALBUILDER_H
#define OHMMS_QMC_GENERIC_MOLECULARORBITALBUILDER_H

#include "OhmmsData/OhmmsElementBase.h"
#include "QMCWaveFunctions/OrbitalBuilderBase.h"

namespace ohmmsqmc {

  class ParticleSet;

  /**
   *@brief Builder class to create Fermi wavefunctions using MO Basis.
   *
   *Any molecular orbital has to use the distance table between the 
   *ions and electrons.
   *The traits of MO Basis is that each basis function is
   *\f$ \phi_{nlm} =  R_{nl} \Re (Y_{lm}) \f$
   *where \f$ Rnl \f$ is any radial grid function and 
   *\f$ \Re (Y_{lm}) \f$ is a real spherical harmonic.
   *Depending upon the input, one can convert the functions on the radial
   *grids or can carry on the calculations using the input functions.
   */
  struct MolecularOrbitalBuilder: public OrbitalBuilderBase {

    typedef map<string,ParticleSet*> PtclPoolType;
    
    MolecularOrbitalBuilder(ParticleSet& p, TrialWaveFunction& psi,
       PtclPoolType& psets):OrbitalBuilderBase(p,psi),ptclPool(psets){ }
  
    bool put(xmlNodePtr cur);
    bool putSpecial(xmlNodePtr cur);
    bool putOpen(const string& fname_in);

    ///need ParticleSetPool
    PtclPoolType& ptclPool;
  };


}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
