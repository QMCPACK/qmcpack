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
/**@file OrbitalBuilderBase.h
 *@brief declaration of the base class for many-body wavefunction.
 */
#ifndef OHMMS_QMC_TRIALORBITALBUILDERBASE_H
#define OHMMS_QMC_TRIALORBITALBUILDERBASE_H 

#include "QMCWaveFunctions/TrialWaveFunction.h"

namespace ohmmsqmc {

  /** An abstract class to build a many-body wavefunction */
  class OrbitalBuilderBase {
    
  protected:

    /// reference to the many-body wavefucntion to which each derived class add a term
    TrialWaveFunction& wfs_ref;

  public:

    typedef TrialWaveFunction::RealType  RealType;    
    typedef TrialWaveFunction::ValueType ValueType;    
    typedef TrialWaveFunction::PosType   PosType;    
    typedef TrialWaveFunction::GradType  GradType;    

    /**
     *\defgroup WavefunctionTags tags for the elements associated with the many-body wavefunctions
     */
    //@{
    /// the element name for a wavefunction
    static string wfs_tag;  
    /// the element name for a parameter
    static string param_tag;  
    /// the element name for a distancetable
    static string dtable_tag;  
    /// the element name for jatrow
    static string jastrow_tag;
    /// the element name for a set of Slater determinants, contains 1..* Slater determinants
    static string detset_tag;
    /// the element name for a Slater determinant, contains 1..* determinants
    static string sd_tag;
    /// the element name for a determinant, may contain (0..*) orbital or parameter element
    static string det_tag;
    /// the element name for single-particle orbital
    static string spo_tag;
    /// the element name for the basis set: basis set contains multiple basis elements
    static string basisset_tag;
    /// the element name for a group of basis functions: basis contains multiple basisfunc elements
    static string basis_tag;
    /// the element name for a basis function
    static string basisfunc_tag;
    //@}

    /// constructor
    OrbitalBuilderBase(TrialWaveFunction& a): wfs_ref(a) { }

    /// process a xml node at cur
    virtual bool put(xmlNodePtr cur) = 0;
  };
  
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/
