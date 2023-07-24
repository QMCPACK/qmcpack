//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//
// File created by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//////////////////////////////////////////////////////////////////////////////////////

#include "QMCWaveFunctions/BsplineFactory/createBsplineReaderT.h"

#include "CPU/SIMD/vmath.hpp"
#include "CPU/e2iphi.h"
#include "QMCWaveFunctions/BsplineFactory/BsplineReaderBaseT.h"
#include "QMCWaveFunctions/BsplineFactory/BsplineSetT.h"
#include "QMCWaveFunctions/BsplineFactory/HybridRepCplxT.h"
#include "QMCWaveFunctions/BsplineFactory/HybridRepRealT.h"
#include "QMCWaveFunctions/BsplineFactory/HybridRepSetReaderT.h"
#include "QMCWaveFunctions/BsplineFactory/SplineC2COMPTargetT.h"
#include "QMCWaveFunctions/BsplineFactory/SplineC2CT.h"
#include "QMCWaveFunctions/BsplineFactory/SplineC2ROMPTargetT.h"
#include "QMCWaveFunctions/BsplineFactory/SplineC2RT.h"
#include "QMCWaveFunctions/BsplineFactory/SplineR2RT.h"
#include "QMCWaveFunctions/BsplineFactory/SplineSetReaderT.h"
#include "QMCWaveFunctions/EinsplineSetBuilderT.h"
#include "QMCWaveFunctions/BsplineFactory/einspline_helper.hpp"
#include "Utilities/ProgressReportEngine.h"
#include <PlatformSelector.hpp>
#include <fftw3.h>

namespace qmcplusplus
{
template<typename T>
struct CreateComplexHelper
{
  static inline std::unique_ptr<BsplineReaderBaseT<T>> createDouble(EinsplineSetBuilderT<T>* e,
                                                                    bool hybrid_rep,
                                                                    const std::string& useGPU)
  {
    using RealType = typename EinsplineSetBuilderT<T>::RealType;
    std::unique_ptr<BsplineReaderBaseT<T>> aReader;

    app_summary() << "    Using real valued spline SPOs with complex double "
                     "precision storage (C2R)."
                  << std::endl;
    if (CPUOMPTargetSelector::selectPlatform(useGPU) == PlatformKind::OMPTARGET)
    {
      app_summary() << "    Running OpenMP offload code path." << std::endl;
      if (hybrid_rep)
      {
        app_summary() << "    Using hybrid orbital representation." << std::endl;
        aReader = std::make_unique<HybridRepSetReaderT<HybridRepCplxT<SplineC2ROMPTargetT<double, T>>>>(e);
      }
      else
        aReader = std::make_unique<SplineSetReaderT<SplineC2ROMPTargetT<double, T>>>(e);
    }
    else
    {
      app_summary() << "    Running on CPU." << std::endl;
      if (hybrid_rep)
      {
        app_summary() << "    Using hybrid orbital representation." << std::endl;
        aReader = std::make_unique<HybridRepSetReaderT<HybridRepCplxT<SplineC2RT<double, T>>>>(e);
      }
      else
        aReader = std::make_unique<SplineSetReaderT<SplineC2RT<double, T>>>(e);
    }

    return aReader;
  }

  static inline std::unique_ptr<BsplineReaderBaseT<T>> createSingle(EinsplineSetBuilderT<T>* e,
                                                                    bool hybrid_rep,
                                                                    const std::string& useGPU)
  {
    using RealType = typename EinsplineSetBuilderT<T>::RealType;
    std::unique_ptr<BsplineReaderBaseT<T>> aReader;

    app_summary() << "    Using real valued spline SPOs with complex single "
                     "precision storage (C2R)."
                  << std::endl;
    if (CPUOMPTargetSelector::selectPlatform(useGPU) == PlatformKind::OMPTARGET)
    {
      app_summary() << "    Running OpenMP offload code path." << std::endl;
      if (hybrid_rep)
      {
        app_summary() << "    Using hybrid orbital representation." << std::endl;
        aReader = std::make_unique<HybridRepSetReaderT<HybridRepCplxT<SplineC2ROMPTargetT<float, T>>>>(e);
      }
      else
        aReader = std::make_unique<SplineSetReaderT<SplineC2ROMPTargetT<float, T>>>(e);
    }
    else
    {
      app_summary() << "    Running on CPU." << std::endl;
      if (hybrid_rep)
      {
        app_summary() << "    Using hybrid orbital representation." << std::endl;
        aReader = std::make_unique<HybridRepSetReaderT<HybridRepCplxT<SplineC2RT<float, T>>>>(e);
      }
      else
        aReader = std::make_unique<SplineSetReaderT<SplineC2RT<float, T>>>(e);
    }

    return aReader;
  }
};

template<typename T>
struct CreateComplexHelper<std::complex<T>>
{
  using ValueType = std::complex<T>;
  using RealType  = typename EinsplineSetBuilderT<ValueType>::RealType;

  static inline std::unique_ptr<BsplineReaderBaseT<ValueType>> createDouble(EinsplineSetBuilderT<ValueType>* e,
                                                                            bool hybrid_rep,
                                                                            const std::string& useGPU)
  {
    std::unique_ptr<BsplineReaderBaseT<ValueType>> aReader;

    app_summary() << "    Using complex valued spline SPOs with complex double "
                     "precision storage (C2C)."
                  << std::endl;
    if (CPUOMPTargetSelector::selectPlatform(useGPU) == PlatformKind::OMPTARGET)
    {
      app_summary() << "    Running OpenMP offload code path." << std::endl;
      if (hybrid_rep)
      {
        app_summary() << "    Using hybrid orbital representation." << std::endl;
        aReader = std::make_unique<HybridRepSetReaderT<HybridRepCplxT<SplineC2COMPTargetT<double, ValueType>>>>(e);
      }
      else
        aReader = std::make_unique<SplineSetReaderT<SplineC2COMPTargetT<double, ValueType>>>(e);
    }
    else
    {
      app_summary() << "    Running on CPU." << std::endl;
      if (hybrid_rep)
      {
        app_summary() << "    Using hybrid orbital representation." << std::endl;
        aReader = std::make_unique<HybridRepSetReaderT<HybridRepCplxT<SplineC2CT<double, ValueType>>>>(e);
      }
      else
        aReader = std::make_unique<SplineSetReaderT<SplineC2CT<double, ValueType>>>(e);
    }

    return aReader;
  }

  static inline std::unique_ptr<BsplineReaderBaseT<ValueType>> createSingle(EinsplineSetBuilderT<ValueType>* e,
                                                                            bool hybrid_rep,
                                                                            const std::string& useGPU)
  {
    std::unique_ptr<BsplineReaderBaseT<ValueType>> aReader;

    app_summary() << "    Using complex valued spline SPOs with complex single "
                     "precision storage (C2C)."
                  << std::endl;
    if (CPUOMPTargetSelector::selectPlatform(useGPU) == PlatformKind::OMPTARGET)
    {
      app_summary() << "    Running OpenMP offload code path." << std::endl;
      if (hybrid_rep)
      {
        app_summary() << "    Using hybrid orbital representation." << std::endl;
        aReader = std::make_unique<HybridRepSetReaderT<HybridRepCplxT<SplineC2COMPTargetT<float, ValueType>>>>(e);
      }
      else
        aReader = std::make_unique<SplineSetReaderT<SplineC2COMPTargetT<float, ValueType>>>(e);
    }
    else
    {
      app_summary() << "    Running on CPU." << std::endl;
      if (hybrid_rep)
      {
        app_summary() << "    Using hybrid orbital representation." << std::endl;
        aReader = std::make_unique<HybridRepSetReaderT<HybridRepCplxT<SplineC2CT<float, ValueType>>>>(e);
      }
      else
        aReader = std::make_unique<SplineSetReaderT<SplineC2CT<float, ValueType>>>(e);
    }

    return aReader;
  }
};

template<typename T>
std::unique_ptr<BsplineReaderBaseT<T>> createBsplineComplexDoubleT(EinsplineSetBuilderT<T>* e,
                                                                   bool hybrid_rep,
                                                                   const std::string& useGPU)
{
  return CreateComplexHelper<T>::createDouble(e, hybrid_rep, useGPU);
}
#ifdef QMC_COMPLEX
template std::unique_ptr<BsplineReaderBaseT<std::complex<float>>> createBsplineComplexDoubleT<std::complex<float>>(
    EinsplineSetBuilderT<std::complex<float>>* e,
    bool hybrid_rep,
    const std::string& useGPU);

template std::unique_ptr<BsplineReaderBaseT<std::complex<double>>> createBsplineComplexDoubleT<std::complex<double>>(
    EinsplineSetBuilderT<std::complex<double>>* e,
    bool hybrid_rep,
    const std::string& useGPU);

#endif
template std::unique_ptr<BsplineReaderBaseT<float>> createBsplineComplexDoubleT<float>(EinsplineSetBuilderT<float>* e,
                                                                                       bool hybrid_rep,
                                                                                       const std::string& useGPU);

template std::unique_ptr<BsplineReaderBaseT<double>> createBsplineComplexDoubleT<double>(
    EinsplineSetBuilderT<double>* e,
    bool hybrid_rep,
    const std::string& useGPU);

template<typename T>
std::unique_ptr<BsplineReaderBaseT<T>> createBsplineComplexSingleT(EinsplineSetBuilderT<T>* e,
                                                                   bool hybrid_rep,
                                                                   const std::string& useGPU)
{
  return CreateComplexHelper<T>::createSingle(e, hybrid_rep, useGPU);
}

#ifdef QMC_COMPLEX
template std::unique_ptr<BsplineReaderBaseT<std::complex<float>>> createBsplineComplexSingleT<std::complex<float>>(
    EinsplineSetBuilderT<std::complex<float>>* e,
    bool hybrid_rep,
    const std::string& useGPU);

template std::unique_ptr<BsplineReaderBaseT<std::complex<double>>> createBsplineComplexSingleT<std::complex<double>>(
    EinsplineSetBuilderT<std::complex<double>>* e,
    bool hybrid_rep,
    const std::string& useGPU);
#endif

template std::unique_ptr<BsplineReaderBaseT<float>> createBsplineComplexSingleT<float>(EinsplineSetBuilderT<float>* e,
                                                                                       bool hybrid_rep,
                                                                                       const std::string& useGPU);

template std::unique_ptr<BsplineReaderBaseT<double>> createBsplineComplexSingleT<double>(
    EinsplineSetBuilderT<double>* e,
    bool hybrid_rep,
    const std::string& useGPU);

template<typename T>
std::unique_ptr<BsplineReaderBaseT<T>> createBsplineRealDoubleT(EinsplineSetBuilderT<T>* e,
                                                                bool hybrid_rep,
                                                                const std::string& useGPU)
{
  app_summary() << "    Using real valued spline SPOs with real double "
                   "precision storage (R2R)."
                << std::endl;
  if (CPUOMPTargetSelector::selectPlatform(useGPU) == PlatformKind::OMPTARGET)
    app_summary() << "OpenMP offload has not been implemented to support "
                     "real valued spline SPOs with real storage!"
                  << std::endl;
  app_summary() << "    Running on CPU." << std::endl;

  std::unique_ptr<BsplineReaderBaseT<T>> aReader;
  if (hybrid_rep)
  {
    app_summary() << "    Using hybrid orbital representation." << std::endl;
    aReader = std::make_unique<HybridRepSetReaderT<HybridRepRealT<SplineR2RT<double, T>>>>(e);
  }
  else
    aReader = std::make_unique<SplineSetReaderT<SplineR2RT<double, T>>>(e);
  return aReader;
}

template std::unique_ptr<BsplineReaderBaseT<float>> createBsplineRealDoubleT<float>(EinsplineSetBuilderT<float>* e,
                                                                                    bool hybrid_rep,
                                                                                    const std::string& useGPU);

template std::unique_ptr<BsplineReaderBaseT<double>> createBsplineRealDoubleT<double>(EinsplineSetBuilderT<double>* e,
                                                                                      bool hybrid_rep,
                                                                                      const std::string& useGPU);

template<typename T>
std::unique_ptr<BsplineReaderBaseT<T>> createBsplineRealSingleT(EinsplineSetBuilderT<T>* e,
                                                                bool hybrid_rep,
                                                                const std::string& useGPU)
{
  app_summary() << "    Using real valued spline SPOs with real single "
                   "precision storage (R2R)."
                << std::endl;
  if (CPUOMPTargetSelector::selectPlatform(useGPU) == PlatformKind::OMPTARGET)
    app_summary() << "OpenMP offload has not been implemented to support "
                     "real valued spline SPOs with real storage!"
                  << std::endl;
  app_summary() << "    Running on CPU." << std::endl;

  std::unique_ptr<BsplineReaderBaseT<T>> aReader;
  if (hybrid_rep)
  {
    app_summary() << "    Using hybrid orbital representation." << std::endl;
    aReader = std::make_unique<HybridRepSetReaderT<HybridRepRealT<SplineR2RT<float, T>>>>(e);
  }
  else
    aReader = std::make_unique<SplineSetReaderT<SplineR2RT<float, T>>>(e);
  return aReader;
}

template std::unique_ptr<BsplineReaderBaseT<float>> createBsplineRealSingleT<float>(EinsplineSetBuilderT<float>* e,
                                                                                    bool hybrid_rep,
                                                                                    const std::string& useGPU);

template std::unique_ptr<BsplineReaderBaseT<double>> createBsplineRealSingleT<double>(EinsplineSetBuilderT<double>* e,
                                                                                      bool hybrid_rep,
                                                                                      const std::string& useGPU);

} // namespace qmcplusplus
