//////////////////////////////////////////////////////////////////
// (c) Copyright 1998-2002 by Jeongnim Kim
//
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
#ifndef OHMMS_RANDOMGENERATOR
#define OHMMS_RANDOMGENERATOR
#include <math.h>

#include <vector>
using std::vector;
#include "ohmms-config.h"

#ifdef HAVE_LIBBOOST
#include "Utilities/BoostRandom.h"
typedef BoostRandom RandomGenerator_t;
#else
#ifdef USE_SPRNG
#include "Utilities/SprngRandom.h"
typedef SprngRandom<0> RandomGenerator_t;
#else
#include "Utilities/RandRandom.h"
typedef RandRandom RandomGenerator_t;
#endif
#endif
extern RandomGenerator_t Random;

#endif

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
