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
#ifndef OHMMS_QMC_TRAITS_H
#define OHMMS_QMC_TRAITS_H

#include "ohmms-config.h"
#include <string>
#include <vector>
#include <map>
#include "OhmmsPETE/TinyVector.h"
#include "OhmmsPETE/Tensor.h"
#include "OhmmsData/OhmmsElementBase.h"
#include "Lattice/Uniform3DGridLayout.h"
#include "ParticleBase/ParticleAttrib.h"
#include "ParticleBase/ParticleBase.h"

namespace ohmmsqmc {

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

  /** traits for QMC variables
   *
   *typedefs for the QMC data types
   */
  struct QMCTraits {
    enum {DIM = OHMMS_DIM};
    typedef OHMMS_INDEXTYPE           IndexType;
    typedef OHMMS_PRECISION           RealType;
    typedef OHMMS_PRECISION           ValueType;
    typedef TinyVector<RealType,DIM>  PosType;
    typedef TinyVector<ValueType,DIM> GradType;
  };

  /** Particle traits to use UniformGridLayout for the ParticleLayout.
   */
  struct PtclOnLatticeTraits {

    //typedef UniformGridLayout<OHMMS_PRECISION,OHMMS_DIM> ParticleLayout_t;
    typedef Uniform3DGridLayout                          ParticleLayout_t;

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
