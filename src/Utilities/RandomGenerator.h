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
#include <cmath>
#include <vector>
using std::vector;
#include "config.h"

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
#include "Utilities/RandRandom.h"
namespace qmcplusplus {
  typedef RandRandom RandomGenerator_t;
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
