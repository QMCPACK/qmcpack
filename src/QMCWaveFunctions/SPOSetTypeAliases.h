//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2018 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_SPOSETTYPEALIASES_H
#define QMCPLUSPLUS_SPOSETTYPEALIASES_H

#include "QMCWaveFunctions/OrbitalSetTraits.h"
#include "OhmmsPETE/OhmmsArray.h"
#include "Particle/ParticleSet.h"
#include "Particle/VirtualParticleSet.h"

namespace qmcplusplus
{

/*! Collects common types for SPOSet consumers and children
 *  When MI and the idea of gaining typdefs and typealiases through
 *  inheritance ambiguous name resolution is soon to follow.
 *  This prevents that and can help make more explicit where
 *  types actually originate.
 */
struct SPOSetTypeAliases
{
  using ValueType = QMCTraits::ValueType;
  typedef OrbitalSetTraits<ValueType>::IndexVector_t IndexVector_t;
  typedef OrbitalSetTraits<ValueType>::ValueVector_t ValueVector_t;
  typedef OrbitalSetTraits<ValueType>::ValueAlignedVector_t ValueAlignedVector_t;
  typedef OrbitalSetTraits<ValueType>::ValueMatrix_t ValueMatrix_t;
  typedef OrbitalSetTraits<ValueType>::GradVector_t  GradVector_t;
  typedef OrbitalSetTraits<ValueType>::GradMatrix_t  GradMatrix_t;
  typedef OrbitalSetTraits<ValueType>::HessVector_t  HessVector_t;
  typedef OrbitalSetTraits<ValueType>::HessMatrix_t  HessMatrix_t;
  typedef OrbitalSetTraits<ValueType>::HessType      HessType;
  typedef Array<HessType,OHMMS_DIM>                  HessArray_t;
  typedef OrbitalSetTraits<ValueType>::GradHessType  GGGType;
  typedef OrbitalSetTraits<ValueType>::GradHessVector_t GGGVector_t;
  typedef OrbitalSetTraits<ValueType>::GradHessMatrix_t GGGMatrix_t;
  typedef OrbitalSetTraits<ValueType>::VGLVector_t      VGLVector_t;
  typedef ParticleSet::Walker_t                      Walker_t;
};

}
#endif
