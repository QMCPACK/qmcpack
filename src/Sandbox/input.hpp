//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by:
//
// File created by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//////////////////////////////////////////////////////////////////////////////////////
// -*- C++ -*-

#ifndef QMCPLUSPLUS_MINIAPPS_INPUT_H
#define QMCPLUSPLUS_MINIAPPS_INPUT_H

#include "ParticleIOUtility.h"
#if defined(USE_NIO)
#include "input/nio.hpp"
#else
#include "input/graphite.hpp"
#endif
#endif
