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
#include "QMCWaveFunctions/EinsplineAdoptor.h"
#include "QMCWaveFunctions/SplineR2RAdoptor.h"
#include <fftw3.h>
#include <QMCWaveFunctions/einspline_helper.hpp>
#include "QMCWaveFunctions/BsplineReaderBase.h"
#include "QMCWaveFunctions/SplineAdoptorReaderP.h"

namespace qmcplusplus
{

  BsplineReaderBase* createBsplineRealDouble(EinsplineSetBuilder* e, int celltype)
  {
    typedef OHMMS_PRECISION RealType;
    BsplineReaderBase* aReader=nullptr;
    if(celltype == SUPERCELL_OPEN)
      aReader= new SplineMixedAdoptorReader<SplineOpenAdoptor<double,RealType,3> >(this);
    else if(celltype == SUPERCELL_SLAB)
      aReader= new SplineMixedAdoptorReader<SplineMixedAdoptor<double,RealType,3> >(this);
    else
      aReader= new SplineAdoptorReader<SplineR2RAdoptor<double,RealType,3> >(this);
    return aReader;
  }
}
