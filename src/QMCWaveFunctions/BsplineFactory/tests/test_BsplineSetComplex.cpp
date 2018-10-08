//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2018 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#include "catch.hpp"

#include "einspline/bspline_base.h"
#include "einspline/bspline_create.h"
#include "einspline/bspline_eval_d.h"
#include "einspline/multi_bspline_create.h"
#include "einspline/multi_bspline_eval_d.h"
#include "einspline/multi_bspline_eval_s.h"

#include "Configuration.h"
#include "QMCWaveFunctions/BsplineFactory/BsplineSet.h"
#include "QMCWaveFunctions/BsplineFactory/SplineC2CAdoptor.h"
#include "Particle/ParticleSet.h"
#include "Batching.h"
#include <iostream>

namespace qmcplusplus
{

TEST_CASE("BsplineSet_Instantiation", "[wavefunction]")
{
  qmcplusplus::BsplineSet<SplineC2CAdoptor<double, QMCTraits::ComplexType>, Batching::SINGLE> SOA_CPU_BSS;
}

}

