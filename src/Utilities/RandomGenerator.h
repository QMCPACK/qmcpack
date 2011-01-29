//////////////////////////////////////////////////////////////////
// (c) Copyright 1998-2002,2003- by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
/** @file RandomGenerator.h
 * @brief Declare a global Random Number Generator
 *
 * Selected among
 * - boost::random
 * - sprng
 * - math::random
 * qmcplusplus::Random() returns a random number [0,1)
 * For OpenMP is enabled, it is important to use thread-safe boost::random. Each
 * thread uses its own random number generator with a distinct seed. This prevents
 * a use of global lock which will slow down the applications significantly.
 */
#ifndef OHMMS_RANDOMGENERATOR
#define OHMMS_RANDOMGENERATOR
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <cmath>
#include <ctime>        

#include <stdint.h>

inline uint32_t make_seed(int i, int n)
{
  return static_cast<uint32_t>(std::time(0))%10474949+(i+1)*n+i;
}

#ifdef HAVE_LIBBOOST

#include "Utilities/BoostRandom.h"
namespace qmcplusplus {
  typedef BoostRandom<OHMMS_PRECISION> RandomGenerator_t;
  extern RandomGenerator_t Random;
}
#else

#ifdef USE_SPRNG
#include "Utilities/SprngRandom.h"
namespace qmcplusplus {
  typedef SprngRandom<0> RandomGenerator_t;
  extern RandomGenerator_t Random;
}
#else
#include "Utilities/SimpleRandom.h"
namespace qmcplusplus {
  typedef SimpleRandom<MTRand> RandomGenerator_t;
  extern RandomGenerator_t Random;
}
#endif
#endif
#endif

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
