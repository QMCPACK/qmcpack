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
#include "CPU/e2iphi.h"
#include "CPU/SIMD/vmath.hpp"
#include "Utilities/ProgressReportEngine.h"
#include "QMCWaveFunctions/EinsplineSetBuilder.h"
#include "QMCWaveFunctions/BsplineFactory/BsplineSet.h"
#include "QMCWaveFunctions/BsplineFactory/SplineC2R.h"
#include "QMCWaveFunctions/BsplineFactory/SplineC2C.h"
#if defined(ENABLE_OFFLOAD)
#include "QMCWaveFunctions/BsplineFactory/SplineC2ROMPTarget.h"
#include "QMCWaveFunctions/BsplineFactory/SplineC2COMPTarget.h"
#endif
#include "QMCWaveFunctions/BsplineFactory/HybridRepCplx.h"
#include <fftw3.h>
#include "QMCWaveFunctions/einspline_helper.hpp"
#include "QMCWaveFunctions/BsplineFactory/BsplineReaderBase.h"
#include "QMCWaveFunctions/BsplineFactory/SplineSetReader.h"
#include "QMCWaveFunctions/BsplineFactory/HybridRepSetReader.h"

namespace qmcplusplus
{
std::unique_ptr<BsplineReaderBase> createBsplineComplexSingle(EinsplineSetBuilder* e,
                                                              bool hybrid_rep,
                                                              const std::string& useGPU)
{
  using RealType = OHMMS_PRECISION;
  std::unique_ptr<BsplineReaderBase> aReader;

#if defined(QMC_COMPLEX)
  app_summary() << "    Using complex valued spline SPOs with complex single precision storage (C2C)." << std::endl;
#if defined(ENABLE_OFFLOAD)
  if (useGPU == "yes")
  {
    if (hybrid_rep)
    {
      app_summary() << "OpenMP offload has not been implemented to support hybrid orbital representation!"
                    << "    Running on CPU." << std::endl;
      app_summary() << "    Using hybrid orbital representation." << std::endl;
      aReader = std::make_unique<HybridRepSetReader<HybridRepCplx<SplineC2C<float>>>>(e);
    }
    else
    {
      app_summary() << "    Running OpenMP offload code path." << std::endl;
      aReader = std::make_unique<SplineSetReader<SplineC2COMPTarget<float>>>(e);
    }
  }
  else
#endif
  {
    app_summary() << "    Running on CPU." << std::endl;
    if (hybrid_rep)
    {
      app_summary() << "    Using hybrid orbital representation." << std::endl;
      aReader = std::make_unique<HybridRepSetReader<HybridRepCplx<SplineC2C<float>>>>(e);
    }
    else
      aReader = std::make_unique<SplineSetReader<SplineC2C<float>>>(e);
  }
#else //QMC_COMPLEX
  app_summary() << "    Using real valued spline SPOs with complex single precision storage (C2R)." << std::endl;
#if defined(ENABLE_OFFLOAD)
  if (useGPU == "yes")
  {
    if (hybrid_rep)
    {
      app_summary() << "OpenMP offload has not been implemented to support hybrid orbital representation!"
                    << "    Running on CPU." << std::endl;
      app_summary() << "    Using hybrid orbital representation." << std::endl;
      aReader = std::make_unique<HybridRepSetReader<HybridRepCplx<SplineC2R<float>>>>(e);
    }
    else
    {
      app_summary() << "    Running OpenMP offload code path." << std::endl;
      aReader = std::make_unique<SplineSetReader<SplineC2ROMPTarget<float>>>(e);
    }
  }
  else
#endif
  {
    app_summary() << "    Running on CPU." << std::endl;
    if (hybrid_rep)
    {
      app_summary() << "    Using hybrid orbital representation." << std::endl;
      aReader = std::make_unique<HybridRepSetReader<HybridRepCplx<SplineC2R<float>>>>(e);
    }
    else
      aReader = std::make_unique<SplineSetReader<SplineC2R<float>>>(e);
  }
#endif
  return aReader;
}
} // namespace qmcplusplus
