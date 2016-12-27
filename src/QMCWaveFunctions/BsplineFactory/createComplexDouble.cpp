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
    
#include "QMCWaveFunctions/BsplineFactory/macro.h"
#include "Numerics/e2iphi.h"
#include "simd/vmath.hpp"
#include "qmc_common.h"
#include <Utilities/ProgressReportEngine.h>
#include "QMCWaveFunctions/EinsplineSetBuilder.h"
#include "QMCWaveFunctions/EinsplineAdoptor.h"
#include "QMCWaveFunctions/SplineC2XAdoptor.h"
#if defined(QMC_ENABLE_SOA_DET)
#include "QMCWaveFunctions/SplineC2RAdoptor.h"
#include "QMCWaveFunctions/SplineC2CAdoptor.h"
#endif
#include <fftw3.h>
#include <QMCWaveFunctions/einspline_helper.hpp>
#include "QMCWaveFunctions/BsplineReaderBase.h"
#include "QMCWaveFunctions/SplineAdoptorReaderP.h"

namespace qmcplusplus
{

  BsplineReaderBase* createBsplineComplexDouble(EinsplineSetBuilder* e, int numOrbs)
  {
    typedef OHMMS_PRECISION RealType;
    BsplineReaderBase* aReader=nullptr;

#if defined(QMC_COMPLEX)

  #if defined(QMC_ENABLE_SOA_DET)
    if(numOrbs>1)
      aReader= new SplineAdoptorReader<SplineC2CSoA<double,RealType> >(e);
    else
  #endif
      aReader= new SplineAdoptorReader<SplineC2CPackedAdoptor<double,RealType,3> >(e);
#else //QMC_COMPLEX

  #if defined(QMC_ENABLE_SOA_DET)
    if(numOrbs>1)
      aReader= new SplineAdoptorReader<SplineC2RSoA<double,RealType> >(e);
    else
  #endif
      aReader= new SplineAdoptorReader<SplineC2RPackedAdoptor<double,RealType,3> >(e);
#endif

    return aReader;
  }
}

