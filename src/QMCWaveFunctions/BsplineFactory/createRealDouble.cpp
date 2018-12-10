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

#include "qmc_common.h"
#include <Utilities/ProgressReportEngine.h>
#include "QMCWaveFunctions/EinsplineSetBuilder.h"
#include "QMCWaveFunctions/BsplineFactory/BsplineSet.h"
#include "QMCWaveFunctions/BsplineFactory/SplineR2RAdoptor.h"
#include "QMCWaveFunctions/BsplineFactory/HybridRealAdoptor.h"
#include <fftw3.h>
#include <QMCWaveFunctions/einspline_helper.hpp>
#include "QMCWaveFunctions/BsplineFactory/BsplineReaderBase.h"
#include "QMCWaveFunctions/BsplineFactory/SplineAdoptorReaderP.h"
#include "QMCWaveFunctions/BsplineFactory/SplineHybridAdoptorReaderP.h"

namespace qmcplusplus
{

  BsplineReaderBase* createBsplineRealDouble(EinsplineSetBuilder* e, bool hybrid_rep)
  {
    BsplineReaderBase* aReader=nullptr;
    if(hybrid_rep)
      aReader= new SplineHybridAdoptorReader<HybridRealSoA<SplineR2RSoA<double,OHMMS_PRECISION> > >(e);
    else
      aReader= new SplineAdoptorReader<SplineR2RSoA<double,OHMMS_PRECISION> >(e);
    return aReader;
  }
}
