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
#include <PlatformSelector.hpp>
#include "SplineSetReader.h"
#include "HybridRepSetReader.h"

namespace qmcplusplus
{
/** create a reader which handles complex (double size real) splines, C2R or C2C case
 *  spline storage and computation precision is ST
 */
template<typename ST>
std::unique_ptr<BsplineReader> createBsplineComplex(EinsplineSetBuilder* e, bool hybrid_rep, const std::string& useGPU)
{
  using RealType = OHMMS_PRECISION;
  std::unique_ptr<BsplineReader> aReader;

#if defined(QMC_COMPLEX)
  app_summary() << "    Using complex valued spline SPOs with complex " << SplineStoragePrecision<ST>::value
                << " precision storage (C2C)." << std::endl;
  if (CPUOMPTargetSelector::selectPlatform(useGPU) == PlatformKind::OMPTARGET)
  {
    app_summary() << "    Running OpenMP offload code path." << std::endl;
    if (hybrid_rep)
    {
      app_summary() << "    Using hybrid orbital representation." << std::endl;
      aReader = std::make_unique<HybridRepSetReader<HybridRepCplx<SplineC2COMPTarget<ST>>>>(e);
    }
    else
      aReader = std::make_unique<SplineSetReader<SplineC2COMPTarget<ST>>>(e);
  }
  else
  {
    app_summary() << "    Running on CPU." << std::endl;
    if (hybrid_rep)
    {
      app_summary() << "    Using hybrid orbital representation." << std::endl;
      aReader = std::make_unique<HybridRepSetReader<HybridRepCplx<SplineC2C<ST>>>>(e);
    }
    else
      aReader = std::make_unique<SplineSetReader<SplineC2C<ST>>>(e);
  }
#else //QMC_COMPLEX
  app_summary() << "    Using real valued spline SPOs with complex " << SplineStoragePrecision<ST>::value
                << " precision storage (C2R)." << std::endl;
  if (CPUOMPTargetSelector::selectPlatform(useGPU) == PlatformKind::OMPTARGET)
  {
    app_summary() << "    Running OpenMP offload code path." << std::endl;
    if (hybrid_rep)
    {
      app_summary() << "    Using hybrid orbital representation." << std::endl;
      aReader = std::make_unique<HybridRepSetReader<HybridRepCplx<SplineC2ROMPTarget<ST>>>>(e);
    }
    else
      aReader = std::make_unique<SplineSetReader<SplineC2ROMPTarget<ST>>>(e);
  }
  else
  {
    app_summary() << "    Running on CPU." << std::endl;
    if (hybrid_rep)
    {
      app_summary() << "    Using hybrid orbital representation." << std::endl;
      aReader = std::make_unique<HybridRepSetReader<HybridRepCplx<SplineC2R<ST>>>>(e);
    }
    else
      aReader = std::make_unique<SplineSetReader<SplineC2R<ST>>>(e);
  }
#endif
  return aReader;
}

std::unique_ptr<BsplineReader> createBsplineComplex(EinsplineSetBuilder* e,
                                                    bool use_single,
                                                    bool hybrid_rep,
                                                    const std::string& useGPU)
{
  if (use_single)
    return createBsplineComplex<float>(e, hybrid_rep, useGPU);
  else
    return createBsplineComplex<double>(e, hybrid_rep, useGPU);
}

} // namespace qmcplusplus
