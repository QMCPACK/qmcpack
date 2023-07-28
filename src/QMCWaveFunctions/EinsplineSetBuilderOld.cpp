//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "QMCWaveFunctions/EinsplineSetBuilder.h"
#include "QMCWaveFunctions/WaveFunctionComponentBuilder.h"
#include "OhmmsData/AttributeSet.h"
#include "Utilities/Timer.h"
#include "Message/Communicate.h"
#include "Message/CommOperators.h"
#include <vector>
#include "ParticleBase/RandomSeqGenerator.h"
#include "CPU/math.hpp"

#include <array>
#include <string_view>

namespace qmcplusplus
{
bool EinsplineSetBuilder::ReadOrbitalInfo(bool skipChecks)
{
  if (!H5File.open(H5FileName, H5F_ACC_RDONLY))
  {
    app_error() << "Could not open HDF5 file \"" << H5FileName << "\" in EinsplineSetBuilder::ReadOrbitalInfo.\n";
    return false;
  }

  // Read format
  std::string format;
  H5File.read(format, "/format");
  H5File.read(Version, "/version");
  app_log() << "  HDF5 orbital file version " << Version[0] << "." << Version[1] << "." << Version[2] << "\n";
  if (format.find("ES") < format.size())
  {
    Format = ESHDF;
    return ReadOrbitalInfo_ESHDF(skipChecks);
  }

  app_error()
      << "EinsplineSetBuilder::ReadOrbitalInfo too old h5 file which is not in ESHDF format! Regenerate the h5 file";
  return false;
}

} // namespace qmcplusplus
