//////////////////////////////////////////////////////////////////
// (c) Copyright 2012-  by Jeongnim Kim and Ken Esler           //
//////////////////////////////////////////////////////////////////
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
