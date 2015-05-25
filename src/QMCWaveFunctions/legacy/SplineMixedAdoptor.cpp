//////////////////////////////////////////////////////////////////
// (c) Copyright 2012-  by Jeongnim Kim and Ken Esler           //
//////////////////////////////////////////////////////////////////
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
