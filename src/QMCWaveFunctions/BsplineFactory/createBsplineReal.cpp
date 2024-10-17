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
#include "SplineSetReader.h"
#include "HybridRepSetReader.h"

namespace qmcplusplus
{
/** create a reader which handles real splines, R2R case
 *  spline storage and computation precision is ST
 */
template<typename ST>
std::unique_ptr<BsplineReader> createBsplineReal(EinsplineSetBuilder* e, bool hybrid_rep)
{
  app_summary() << "    Using real valued spline SPOs with real " << SplineStoragePrecision<ST>::value
                << " precision storage (R2R)." << std::endl;

  std::unique_ptr<BsplineReader> aReader;
  if (hybrid_rep)
  {
    app_summary() << "    Using hybrid orbital representation." << std::endl;
    aReader = std::make_unique<HybridRepSetReader<HybridRepReal<SplineR2R<ST>>>>(e);
  }
  else
    aReader = std::make_unique<SplineSetReader<SplineR2R<ST>>>(e);
  return aReader;
}

std::unique_ptr<BsplineReader> createBsplineReal(EinsplineSetBuilder* e,
                                                 bool use_single,
                                                 bool hybrid_rep)
{
  if (use_single)
    return createBsplineReal<float>(e, hybrid_rep);
  else
    return createBsplineReal<double>(e, hybrid_rep);
}

} // namespace qmcplusplus
