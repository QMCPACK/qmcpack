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
    
    
/** @file SplineC2XAdoptor.cpp
 * instantiate C2X adoptor and its reader.
 */
#include <Configuration.h>
#include <Numerics/e2iphi.h>
#include <simd/vmath.hpp>
#include <Utilities/ProgressReportEngine.h>
#include <fftw3.h>
#include <Message/CommOperators.h>
#include <QMCWaveFunctions/EinsplineSetBuilder.h>
#include <QMCWaveFunctions/EinsplineAdoptor.h>
#include <spline/einspline_util.hpp>
#include <QMCWaveFunctions/einspline_helper.hpp>
#include <QMCWaveFunctions/SplineC2XAdoptor.h>
#include <QMCWaveFunctions/SplineAdoptorReader.h>
namespace qmcplusplus
{
#if defined(QMC_COMPLEX)
template class SplineC2CPackedAdoptor<float,double,3>;
template class SplineAdoptorReader<SplineC2CPackedAdoptor<float,double,3> >;
#else
template class SplineC2RPackedAdoptor<float,double,3>;
template class SplineAdoptorReader<SplineC2RPackedAdoptor<float,double,3> >;
#endif
}
