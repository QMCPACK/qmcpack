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
#ifndef QMCPLUSPLUS_TRAITS_H
#define QMCPLUSPLUS_TRAITS_H

#include "ohmms-config.h"
#include <string>
#include <vector>
#include <map>
#include <complex>
#include "OhmmsPETE/TinyVector.h"
#include "OhmmsPETE/Tensor.h"
#include "OhmmsData/OhmmsElementBase.h"
#include "Lattice/Uniform3DGridLayout.h"
#include "ParticleBase/ParticleAttrib.h"
#include "ParticleBase/ParticleBase.h"
#include "Utilities/OhmmsInfo.h"
#include "Message/Communicate.h"

namespace qmcplusplus {

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
    typedef OHMMS_INDEXTYPE                IndexType;
    typedef OHMMS_PRECISION                RealType;
#if defined(QMC_COMPLEX)
    typedef std::complex<OHMMS_PRECISION>  ValueType;
#else
    typedef OHMMS_PRECISION                ValueType;
#endif
    typedef std::complex<OHMMS_PRECISION>  ComplexType;
    typedef TinyVector<RealType,DIM>       PosType;
    typedef TinyVector<ValueType,DIM>      GradType;
    typedef Tensor<RealType,DIM>           TensorType;

  };

  /** Particle traits to use UniformGridLayout for the ParticleLayout.
   */
  struct PtclOnLatticeTraits{

    //typedef UniformGridLayout<OHMMS_PRECISION,OHMMS_DIM> ParticleLayout_t;
    typedef Uniform3DGridLayout                          ParticleLayout_t;

    typedef int                                          Index_t;
    typedef OHMMS_PRECISION                              Scalar_t;
    typedef std::complex<Scalar_t>                       Complex_t;

    typedef ParticleLayout_t::SingleParticleIndex_t      SingleParticleIndex_t;
    typedef ParticleLayout_t::SingleParticlePos_t        SingleParticlePos_t;
    typedef ParticleLayout_t::Tensor_t                   Tensor_t;

    typedef ParticleAttrib<Index_t>                      ParticleIndex_t;
    typedef ParticleAttrib<Scalar_t>                     ParticleScalar_t;
    typedef ParticleAttrib<SingleParticlePos_t>          ParticlePos_t;
    typedef ParticleAttrib<Tensor_t>                     ParticleTensor_t;

#if defined(QMC_COMPLEX)
    typedef ParticleAttrib<TinyVector<Complex_t,OHMMS_DIM> > ParticleGradient_t;
    typedef ParticleAttrib<Complex_t>                      ParticleLaplacian_t;
#else
    typedef ParticleAttrib<SingleParticlePos_t>            ParticleGradient_t;
    typedef ParticleAttrib<Scalar_t>                       ParticleLaplacian_t;
#endif
  };

  inline std::ostream& app_log() { 
    return  OhmmsInfo::Log->getStream(); 
  }

  inline std::ostream& app_error(){ 
    OhmmsInfo::Log->getStream() << "ERROR ";
    return OhmmsInfo::Error->getStream();
  }
  inline std::ostream& app_warning(){ 
    OhmmsInfo::Log->getStream() << "WARNING ";
    return OhmmsInfo::Warn->getStream();
  }
  inline std::ostream& app_debug(){ 
    OhmmsInfo::Log->getStream() << "DEBUG ";
    return OhmmsInfo::Debug->getStream();
  }
}

#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
