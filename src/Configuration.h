//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_TRAITS_H
#define QMCPLUSPLUS_TRAITS_H

#include <config.h>
#include <string>
#include <vector>
#include <map>
#include "type_traits/QMCTypes.h"
#include "OhmmsData/OhmmsElementBase.h"
#include "OhmmsData/RecordProperty.h"
#include "Particle/Lattice/CrystalLattice.h"
#include "Particle/ParticleBase/ParticleAttrib.h"
#include "Platforms/Host/OutputManager.h"
#include "Message/Communicate.h"

//define empty DEBUG_MEMORY
#define DEBUG_MEMORY(msg)
//uncomment this out to trace the call tree of destructors
//#define DEBUG_MEMORY(msg) std::cerr << "<<<< " << msg << std::endl;

#if defined(DEBUG_PSIBUFFER_ON)
#define DEBUG_PSIBUFFER(who, msg)                              \
  std::cerr << "PSIBUFFER " << who << " " << msg << std::endl; \
  std::cerr.flush();
#else
#define DEBUG_PSIBUFFER(who, msg)
#endif

namespace qmcplusplus
{
/** traits for QMC variables
 *
 *typedefs for the QMC data types
 */
struct QMCTraits
{
  enum
  {
    DIM = OHMMS_DIM
  };
  using QTBase = QMCTypes<OHMMS_PRECISION, DIM>;
  using QTFull = QMCTypes<OHMMS_PRECISION_FULL, DIM>;
  typedef QTBase::RealType RealType;
  typedef QTBase::ComplexType ComplexType;
  typedef QTBase::ValueType ValueType;
  typedef QTBase::PosType PosType;
  typedef QTBase::GradType GradType;
  typedef QTBase::TensorType TensorType;
  ///define other types
  typedef OHMMS_INDEXTYPE IndexType;
  typedef QTFull::RealType FullPrecRealType;
  typedef QTFull::ValueType FullPrecValueType;
  ///define PropertyList_t
  typedef RecordNamedProperty<FullPrecRealType> PropertySetType;

  // Type for particle group index pairs
  using PtclGrpIndexes = std::vector<std::pair<int,int>>;
};

/** Particle traits to use UniformGridLayout for the ParticleLayout.
 */
struct PtclOnLatticeTraits
{
  using ParticleLayout = CrystalLattice<OHMMS_PRECISION, OHMMS_DIM>;
  using QTFull = QMCTraits::QTFull;

  typedef int Index_t;
  typedef QTFull::RealType Scalar_t;
  typedef QTFull::ComplexType Complex_t;

  typedef ParticleLayout::SingleParticleIndex SingleParticleIndex;
  typedef ParticleLayout::SingleParticlePos SingleParticlePos;
  typedef ParticleLayout::Tensor_t Tensor_t;

  typedef ParticleAttrib<Index_t> ParticleIndex;
  typedef ParticleAttrib<Scalar_t> ParticleScalar;
  typedef ParticleAttrib<SingleParticlePos> ParticlePos;
  typedef ParticleAttrib<Tensor_t> ParticleTensor;

  typedef ParticleAttrib<QTFull::GradType> ParticleGradient;
  typedef ParticleAttrib<QTFull::ValueType> ParticleLaplacian;
  typedef QTFull::ValueType SingleParticleValue;
};


// For unit tests
//  Check if we are compiling with Catch defined.  Could use other symbols if needed.
#ifdef TEST_CASE
#ifdef QMC_COMPLEX
  using ValueApprox = Catch::Detail::ComplexApprox;
#else
  using ValueApprox = Catch::Detail::Approx;
#endif
#endif

} // namespace qmcplusplus

#endif
