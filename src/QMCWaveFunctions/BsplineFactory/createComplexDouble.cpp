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


#include "QMCWaveFunctions/BsplineFactory/createBsplineReader.h"
#include "Numerics/e2iphi.h"
#include "CPU/SIMD/vmath.hpp"
#include <Utilities/ProgressReportEngine.h>
#include "QMCWaveFunctions/EinsplineSetBuilder.h"
#include "QMCWaveFunctions/BsplineFactory/BsplineSet.h"
#include "QMCWaveFunctions/BsplineFactory/SplineC2R.h"
#include "QMCWaveFunctions/BsplineFactory/SplineC2C.h"
#if defined(ENABLE_OFFLOAD)
#include "QMCWaveFunctions/BsplineFactory/SplineC2ROMP.h"
#endif
#include "QMCWaveFunctions/BsplineFactory/HybridRepCplx.h"
#include <fftw3.h>
#include <QMCWaveFunctions/einspline_helper.hpp>
#include "QMCWaveFunctions/BsplineFactory/BsplineReaderBase.h"
#include "QMCWaveFunctions/BsplineFactory/SplineSetReader.h"
#include "QMCWaveFunctions/BsplineFactory/HybridRepSetReader.h"

namespace qmcplusplus
{
BsplineReaderBase* createBsplineComplexDouble(EinsplineSetBuilder* e, bool hybrid_rep, const std::string& useGPU)
{
  typedef OHMMS_PRECISION RealType;
  BsplineReaderBase* aReader = nullptr;

#if defined(QMC_COMPLEX)
  if (hybrid_rep)
    aReader = new HybridRepSetReader<HybridRepCplx<SplineC2C<double>>>(e);
  else
    aReader = new SplineSetReader<SplineC2C<double>>(e);
#else //QMC_COMPLEX
#if defined(ENABLE_OFFLOAD)
  if (useGPU == "yes")
  {
    if (hybrid_rep)
    {
      APP_ABORT("OpenMP offload has not been enabled with hybrid orbital representation!");
    }
    else
      aReader = new SplineSetReader<SplineC2ROMP<double>>(e);
  }
  else
#endif
  {
    if (hybrid_rep)
      aReader = new HybridRepSetReader<HybridRepCplx<SplineC2R<double>>>(e);
    else
      aReader = new SplineSetReader<SplineC2R<double>>(e);
  }
#endif

  return aReader;
}
} // namespace qmcplusplus
