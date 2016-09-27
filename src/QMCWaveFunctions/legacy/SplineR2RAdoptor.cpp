//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
/** @file SplineR2RAdoptor.cpp
 * instantiate SplineR2RAdoptor and it reader
 */
#include <Configuration.h>
#include <Numerics/e2iphi.h>
#include <simd/vmath.hpp>
#include <fftw3.h>
#include <Utilities/ProgressReportEngine.h>
#include <Message/CommOperators.h>
#include <QMCWaveFunctions/EinsplineSetBuilder.h>
#include <QMCWaveFunctions/einspline_helper.hpp>
#include <QMCWaveFunctions/EinsplineAdoptor.h>
#include <QMCWaveFunctions/SplineR2RAdoptor.h>
#include <QMCWaveFunctions/SplineAdoptorReader.h>
namespace qmcplusplus
{
template class SplineR2RAdoptor<float,double,3>;
template class SplineAdoptorReader<SplineR2RAdoptor<float,double,3> >;
}
