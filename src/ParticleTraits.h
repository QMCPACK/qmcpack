//////////////////////////////////////////////////////////////////
// (c) Copyright 2003- by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim
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
//   Department of Physics, Ohio State University
//   Ohio Supercomputer Center
//////////////////////////////////////////////////////////////////
// -*- C++ -*-

#ifndef OHMMS_PARTICLE_TRAITS_COLLECTION_H
#define OHMMS_PARTICLE_TRAITS_COLLECTION_H

#ifdef HAVE_CONFIG_H
#include "ohmms-config.h"
#endif
#include <string>
#include <vector>
#include <map>
#include "OhmmsPETE/TinyVector.h"
#include "OhmmsPETE/Tensor.h"
#include "Lattice/UniformGridLayout.h"
#include "ParticleBase/ParticleAttrib.h"
#include "ParticleBase/ParticleBase.h"

/**@file ParticleTraits.h
 *@brief Declaration of the trait classe for ParticleBase
 */

namespace OHMMS {

  /** traits for the common particle attributes
   *
   *This is an alternative to the global typedefs.
   */
  struct PtclAttribTraits {

    typedef int                                                     Index_t;
    typedef ParticleAttrib<Index_t>                                 ParticleIndex_t;
    typedef ParticleAttrib<OHMMS_PRECISION>                         ParticleScalar_t;
    typedef ParticleAttrib<TinyVector<OHMMS_PRECISION, OHMMS_DIM> > ParticlePos_t;
    typedef ParticleAttrib<Tensor<OHMMS_PRECISION, OHMMS_DIM> >     ParticleTensor_t;

  };

  /** Particle traits to use UniformGridLayout for the ParticleLayout.
   */
  struct PtclOnLatticeTraits {

    typedef UniformGridLayout<OHMMS_PRECISION,OHMMS_DIM> ParticleLayout_t;

    typedef int                                          Index_t;
    typedef ParticleLayout_t::Scalar_t                   Scalar_t;
    typedef ParticleLayout_t::SingleParticleIndex_t      SingleParticleIndex_t;
    typedef ParticleLayout_t::SingleParticlePos_t        SingleParticlePos_t;
    typedef ParticleLayout_t::Tensor_t                   Tensor_t;

    typedef ParticleAttrib<Index_t>                      ParticleIndex_t;
    typedef ParticleAttrib<Scalar_t>                     ParticleScalar_t;
    typedef ParticleAttrib<SingleParticlePos_t>          ParticlePos_t;
    typedef ParticleAttrib<Tensor_t>                     ParticleTensor_t;


  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/
