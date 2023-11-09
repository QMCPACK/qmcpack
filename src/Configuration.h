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
  using QTBase      = QMCTypes<OHMMS_PRECISION, DIM>;
  using QTFull      = QMCTypes<OHMMS_PRECISION_FULL, DIM>;
  using RealType    = QTBase::RealType;
  using ComplexType = QTBase::ComplexType;
  using ValueType   = QTBase::ValueType;
  using PosType     = QTBase::PosType;
  using GradType    = QTBase::GradType;
  using TensorType  = QTBase::TensorType;
  ///define other types
  using IndexType         = OHMMS_INDEXTYPE;
  using FullPrecRealType  = QTFull::RealType;
  using FullPrecValueType = QTFull::ValueType;
  ///define PropertyList_t
  using PropertySetType = RecordNamedProperty<FullPrecRealType>;

  // Type for particle group index pairs
  using PtclGrpIndexes = std::vector<std::pair<int, int>>;
};

/** Particle traits to use UniformGridLayout for the ParticleLayout.
 */
struct PtclOnLatticeTraits
{
  using ParticleLayout = CrystalLattice<OHMMS_PRECISION, OHMMS_DIM>;
  using QTFull         = QMCTraits::QTFull;

  using Index_t   = int;
  using Scalar_t  = QTFull::RealType;
  using Complex_t = QTFull::ComplexType;

  using SingleParticleIndex = ParticleLayout::SingleParticleIndex;
  using SingleParticlePos   = ParticleLayout::SingleParticlePos;
  using Tensor_t            = ParticleLayout::Tensor_t;

  using ParticleIndex  = ParticleAttrib<Index_t>;
  using ParticleScalar = ParticleAttrib<Scalar_t>;
  using ParticlePos    = ParticleAttrib<SingleParticlePos>;
  using ParticleTensor = ParticleAttrib<Tensor_t>;

  using ParticleGradient    = ParticleAttrib<QTFull::GradType>;
  using ParticleLaplacian   = ParticleAttrib<QTFull::ValueType>;
  using SingleParticleValue = QTFull::ValueType;
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
