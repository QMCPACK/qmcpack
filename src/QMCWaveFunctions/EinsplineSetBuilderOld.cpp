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
    app_error() << "Could not open HDF5 file \"" << H5FileName
                << "\" in EinsplineSetBuilder::ReadOrbitalInfo.  Aborting.\n";
    APP_ABORT("EinsplineSetBuilder::ReadOrbitalInfo");
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
  //////////////////////////////////////////////////
  // Read basic parameters from the orbital file. //
  //////////////////////////////////////////////////
  // Check the version
  if (Version[0] == 0 && Version[1] == 11)
  {
    parameterGroup   = "/parameters_0";
    ionsGroup        = "/ions_2";
    eigenstatesGroup = "/eigenstates_3";
  }
  else if (Version[0] == 0 && Version[1] == 20)
  {
    parameterGroup   = "/parameters";
    ionsGroup        = "/ions";
    eigenstatesGroup = "/eigenstates";
  }
  else
  {
    std::ostringstream o;
    o << "Unknown HDF5 orbital file version " << Version[0] << "." << Version[1] << "." << Version[2] << "\n";
    APP_ABORT(o.str());
  }
  H5File.read(Lattice, parameterGroup + "/lattice");
  H5File.read(RecipLattice, parameterGroup + "/reciprocal_lattice");
  SuperLattice = dot(TileMatrix, Lattice);
  std::array<char, 1000> buff;
  int length = std::snprintf(buff.data(), buff.size(),
                             "  Lattice = \n    [ %8.5f %8.5f %8.5f\n"
                             "      %8.5f %8.5f %8.5f\n"
                             "      %8.5f %8.5f %8.5f ]\n",
                             Lattice(0, 0), Lattice(0, 1), Lattice(0, 2), Lattice(1, 0), Lattice(1, 1), Lattice(1, 2),
                             Lattice(2, 0), Lattice(2, 1), Lattice(2, 2));
  if (length < 0)
    throw std::runtime_error("Error generating Lattice string");
  app_log() << std::string_view(buff.data(), length);
  length =
      std::snprintf(buff.data(), buff.size(),
                    "  SuperLattice = \n    [ %13.12f %13.12f %13.12f\n"
                    "      %13.12f %13.12f %13.12f\n"
                    "      %13.12f %13.12f %13.12f ]\n",
                    SuperLattice(0, 0), SuperLattice(0, 1), SuperLattice(0, 2), SuperLattice(1, 0), SuperLattice(1, 1),
                    SuperLattice(1, 2), SuperLattice(2, 0), SuperLattice(2, 1), SuperLattice(2, 2));
  if (length < 0)
    throw std::runtime_error("Error generating SuperLattice string");
  if (!CheckLattice())
    APP_ABORT("CheckLattice failed");
  app_log() << std::string_view(buff.data(), length);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      LatticeInv(i, j) = RecipLattice(i, j) / (2.0 * M_PI);
  H5File.read(NumBands, parameterGroup + "/num_bands");
  H5File.read(NumElectrons, parameterGroup + "/num_electrons");
  H5File.read(NumSpins, parameterGroup + "/num_spins");
  H5File.read(NumTwists, parameterGroup + "/num_twists");
  app_log() << "bands=" << NumBands << ", elecs=" << NumElectrons << ", spins=" << NumSpins << ", twists=" << NumTwists
            << std::endl;

  //////////////////////////////////
  // Read ion types and locations //
  //////////////////////////////////
  H5File.read(IonTypes, ionsGroup + "/atom_types");
  H5File.read(IonPos, ionsGroup + "/pos");
  ///////////////////////////
  // Read the twist angles //
  ///////////////////////////
  primcell_kpoints.resize(NumTwists);
  for (int ti = 0; ti < NumTwists; ti++)
  {
    std::ostringstream path;
    if ((Version[0] == 0 && Version[1] == 11) || NumTwists > 1)
      path << eigenstatesGroup << "/twist_" << ti << "/twist_angle";
    else
      path << eigenstatesGroup << "/twist/twist_angle";
    TinyVector<double, OHMMS_DIM> primcell_kpoints_DP;
    H5File.read(primcell_kpoints_DP, path.str());
    primcell_kpoints[ti] = primcell_kpoints_DP;
    int length      = std::snprintf(buff.data(), buff.size(), "  Found twist angle (%6.3f, %6.3f, %6.3f)\n",
                               primcell_kpoints[ti][0], primcell_kpoints[ti][1], primcell_kpoints[ti][2]);
    if (length < 0)
      throw std::runtime_error("Error converting twist angle to string");
    app_log() << std::string_view(buff.data(), length);
  }
  //////////////////////////////////////////////////////////
  // If the density has not been set in TargetPtcl, and   //
  // the density is available, read it in and save it     //
  // in TargetPtcl.                                       //
  //////////////////////////////////////////////////////////
  if (TargetPtcl.Density_G.empty())
  {
    Array<double, OHMMS_DIM> Density_r_DP;
    H5File.read(TargetPtcl.DensityReducedGvecs, "/density/reduced_gvecs");
    H5File.read(Density_r_DP, "/density/rho_r");
    TargetPtcl.Density_r = Density_r_DP;
    int numG             = TargetPtcl.DensityReducedGvecs.size();
    // Convert primitive G-vectors to supercell G-vectors
    for (int iG = 0; iG < numG; iG++)
      TargetPtcl.DensityReducedGvecs[iG] = dot(TileMatrix, TargetPtcl.DensityReducedGvecs[iG]);
    app_log() << "  Read " << numG << " density G-vectors.\n";
    if (TargetPtcl.DensityReducedGvecs.size())
    {
      app_log() << "  EinsplineSetBuilder found density in the HDF5 file.\n";
      std::vector<std::complex<double>> Density_G_DP;
      H5File.read(Density_G_DP, "/density/rho_G");
      TargetPtcl.Density_G.assign(Density_G_DP.begin(), Density_G_DP.end());
      if (!TargetPtcl.Density_G.size())
      {
        app_error() << "  Density reduced G-vectors defined, but not the"
                    << " density.\n";
        abort();
      }
    }
  }
  return true;
}

} // namespace qmcplusplus
