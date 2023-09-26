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

namespace qmcplusplus
{
bool EinsplineSetBuilder::ReadOrbitalInfo(bool skipChecks)
{
  if (!H5File.open(H5FileName, H5F_ACC_RDONLY))
  {
    app_error() << "Could not open HDF5 file \"" << H5FileName << "\" in EinsplineSetBuilder::ReadOrbitalInfo.\n";
    return false;
  }

  try
  {
    // Read format
    std::string format;
    H5File.read(format, "/format");
    if (format.find("ES") == std::string::npos)
      throw std::runtime_error("Format string input \"" + format + "\" doesn't contain \"ES\" keyword.");
    Format = ESHDF;
    H5File.read(Version, "/version");
    app_log() << "  HDF5 orbital file version " << Version[0] << "." << Version[1] << "." << Version[2] << std::endl;
  }
  catch (const std::exception& e)
  {
    app_error() << e.what() << std::endl
                << "EinsplineSetBuilder::ReadOrbitalInfo too old h5 file which is not in ESHDF format!" << std::endl;
    return false;
  }

  return ReadOrbitalInfo_ESHDF(skipChecks);
}

} // namespace qmcplusplus
