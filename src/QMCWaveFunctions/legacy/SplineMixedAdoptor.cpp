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
    
    
/** @file SplineMixedAdoptor.cpp
 * instantiate SplineMixedAdoptor and  reader
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
#include <QMCWaveFunctions/SplineMixedAdoptor.h>
#include <QMCWaveFunctions/SplineMixedAdoptorReader.h>
namespace qmcplusplus
{
template class SplineMixedAdoptor<float,double,3>;
template class SplineMixedAdoptorReader<SplineMixedAdoptor<float,double,3> >;

template class SplineOpenAdoptor<float,double,3>;
template class SplineMixedAdoptorReader<SplineOpenAdoptor<float,double,3> >;
}
