//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "EinsplineSetBuilder.h"
#include "DistanceTable.h"
#include "OhmmsData/AttributeSet.h"
#include "Utilities/Timer.h"
#include "Message/Communicate.h"
#include "Message/CommOperators.h"
#include "ParticleBase/RandomSeqGenerator.h"
#include "Utilities/qmc_common.h"

#include <array>
#include <string_view>

namespace qmcplusplus
{
bool sortByIndex(BandInfo leftB, BandInfo rightB)
{
  if (leftB.BandIndex == rightB.BandIndex)
  {
    if ((leftB.Energy < rightB.Energy + 1e-6) && (leftB.Energy > rightB.Energy - 1e-6))
      return leftB.TwistIndex < rightB.TwistIndex;
    else
      return leftB.Energy < rightB.Energy;
  }
  else
    return (leftB.BandIndex < rightB.BandIndex);
};

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
                << "EinsplineSetBuilder::ReadOrbitalInfo h5 file format is too old or it is not a bspline orbital file!"
                << std::endl;
    return false;
  }

  return ReadOrbitalInfo_ESHDF(skipChecks);
}

bool EinsplineSetBuilder::ReadOrbitalInfo_ESHDF(bool skipChecks)
{
  app_log() << "  Reading orbital file in ESHDF format.\n";
  H5File.read(Version, "/version");
  app_log() << "  ESHDF orbital file version " << Version[0] << "." << Version[1] << "." << Version[2] << std::endl;
  H5File.read(Lattice, "/supercell/primitive_vectors");
  RecipLattice = 2.0 * M_PI * inverse(Lattice);
  SuperLattice = dot(TileMatrix, Lattice);
  std::array<char, 1000> buff;
  int length = std::snprintf(buff.data(), buff.size(),
                             "  Lattice = \n    [ %9.6f %9.6f %9.6f\n"
                             "      %9.6f %9.6f %9.6f\n"
                             "      %9.6f %9.6f %9.6f ]\n",
                             Lattice(0, 0), Lattice(0, 1), Lattice(0, 2), Lattice(1, 0), Lattice(1, 1), Lattice(1, 2),
                             Lattice(2, 0), Lattice(2, 1), Lattice(2, 2));
  if (length < 0)
    throw std::runtime_error("Error converting lattice to a string");
  app_log() << std::string_view(buff.data(), length);
  length =
      std::snprintf(buff.data(), buff.size(),
                    "  SuperLattice = \n    [ %9.6f %9.6f %9.6f\n"
                    "      %9.6f %9.6f %9.6f\n"
                    "      %9.6f %9.6f %9.6f ]\n",
                    SuperLattice(0, 0), SuperLattice(0, 1), SuperLattice(0, 2), SuperLattice(1, 0), SuperLattice(1, 1),
                    SuperLattice(1, 2), SuperLattice(2, 0), SuperLattice(2, 1), SuperLattice(2, 2));
  if (length < 0)
    throw std::runtime_error("Error converting SuperLattice to a string");
  app_log() << std::string_view(buff.data(), length) << std::endl;
  if (!CheckLattice())
    throw std::runtime_error("CheckLattice failed");
  PrimCell.set(Lattice);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      LatticeInv(i, j) = RecipLattice(i, j) / (2.0 * M_PI);
  int have_dpsi = false;
  NumTwists = NumSpins = NumBands = 0;
  NumElectrons                    = TargetPtcl.getTotalNum();
  H5File.read(NumBands, "/electrons/kpoint_0/spin_0/number_of_states");
  H5File.readEntry(NumSpins, "/electrons/number_of_spins");
  H5File.read(NumTwists, "/electrons/number_of_kpoints");
  H5File.readEntry(have_dpsi, "/electrons/have_dpsi");
  HaveOrbDerivs = have_dpsi;
  app_log() << "bands=" << NumBands << ", elecs=" << NumElectrons << ", spins=" << NumSpins << ", twists=" << NumTwists
            << std::endl;
  //////////////////////////////////
  // Read ion types and locations //
  //////////////////////////////////
  Vector<int> species_ids;
  H5File.read(species_ids, "/atoms/species_ids");
  int num_species;
  H5File.read(num_species, "/atoms/number_of_species");
  std::vector<int> atomic_numbers(num_species);
  for (int isp = 0; isp < num_species; isp++)
  {
    std::ostringstream name;
    name << "/atoms/species_" << isp << "/atomic_number";
    H5File.readEntry(atomic_numbers[isp], name.str());
  }
  IonTypes.resize(species_ids.size());
  for (int i = 0; i < species_ids.size(); i++)
    IonTypes[i] = atomic_numbers[species_ids[i]];
  H5File.read(IonPos, "/atoms/positions");
  for (int i = 0; i < IonTypes.size(); i++)
    app_log() << "Atom type(" << i << ") = " << IonTypes[i] << std::endl;
  /////////////////////////////////////
  // Read atom orbital info from xml //
  /////////////////////////////////////
  // construct Super2Prim mapping.
  if (Super2Prim.size() == 0)
  {
    //SourcePtcl->convert2Cart(SourcePtcl->R);
    Super2Prim.resize(SourcePtcl->R.size(), -1);
    std::vector<int> prim_atom_counts;
    prim_atom_counts.resize(IonPos.size(), 0);
    for (int i = 0; i < SourcePtcl->R.size(); i++)
    {
      PosType ref = PrimCell.toUnit_floor(SourcePtcl->R[i]);
      for (int j = 0; j < IonPos.size(); j++)
      {
        PosType dr = PrimCell.toUnit_floor(IonPos[j]) - ref;
        for (int k = 0; k < OHMMS_DIM; k++)
          dr[k] -= round(dr[k]);
        if (dot(dr, dr) < MatchingTol)
        {
          if (Super2Prim[i] < 0)
          {
            Super2Prim[i] = j;
            prim_atom_counts[j]++;
          }
          else
          {
            app_error() << "Supercell ion " << i << " at " << SourcePtcl->R[j]
                        << " was found twice in the primitive cell as ion " << Super2Prim[i] << " and " << j
                        << std::endl;
            if (!skipChecks)
              abort();
          }
        }
      }
      if (Super2Prim[i] < 0)
      {
        app_error() << "Supercell ion " << i << " not found in the primitive cell" << std::endl;
        if (!skipChecks)
          abort();
      }
      else
      {
        //app_log() << "Supercell ion " << i << " mapped to primitive cell ion " << Super2Prim[i] << std::endl;
      }
    }
    const int tiling_size = std::abs(det(TileMatrix));
    for (int i = 0; i < IonPos.size(); i++)
      if (prim_atom_counts[i] != tiling_size)
      {
        app_error() << "Primitive cell ion " << i << " was found only " << prim_atom_counts[i]
                    << " times in the supercell rather than " << tiling_size << std::endl;
        if (!skipChecks)
          abort();
      }
    // construct AtomicCentersInfo
    AtomicCentersInfo.resize(IonPos.size());
    for (int i = 0; i < IonPos.size(); i++)
      AtomicCentersInfo.ion_pos[i] = IonPos[i];
    const auto& source_species = SourcePtcl->getSpeciesSet();
    int Zind                   = source_species.findAttribute("atomicnumber");
    const int table_id         = SourcePtcl->addTable(*SourcePtcl);
    const auto& ii_table       = SourcePtcl->getDistTable(table_id);
    SourcePtcl->update(true);
    for (int i = 0; i < IonPos.size(); i++)
    {
      AtomicCentersInfo.non_overlapping_radius[i] = std::numeric_limits<RealType>::max();
      //should only call get_first_neighbor to set non_overlapping_radius if there are more than one atom  in the cell
      if (Super2Prim.size() == 1)
        continue;
      for (int j = 0; j < Super2Prim.size(); j++)
        if (Super2Prim[j] == i)
        {
          // set GroupID for each ion in primitive cell
          if ((Zind < 0) || (source_species(Zind, SourcePtcl->GroupID[j]) == IonTypes[i]))
            AtomicCentersInfo.GroupID[i] = SourcePtcl->GroupID[j];
          else
          {
            app_error() << "Primitive cell ion " << i << " vs supercell ion " << j
                        << " atomic number not matching: " << IonTypes[i] << " vs "
                        << source_species(Zind, SourcePtcl->GroupID[j]) << std::endl;
            if (!skipChecks)
              abort();
          }
          // set non_overlapping_radius for each ion in primitive cell
          RealType r(0);
          PosType dr;
          ii_table.get_first_neighbor(j, r, dr, false);
          if (r < 1e-3)
            APP_ABORT("EinsplineSetBuilder::ReadOrbitalInfo_ESHDF too close ions <1e-3 bohr!");
          AtomicCentersInfo.non_overlapping_radius[i] = 0.5 * r;
          break;
        }
    }

    // load cutoff_radius, spline_radius, spline_npoints, lmax if exists.
    const int inner_cutoff_ind   = source_species.findAttribute("inner_cutoff");
    const int cutoff_radius_ind  = source_species.findAttribute("cutoff_radius");
    const int spline_radius_ind  = source_species.findAttribute("spline_radius");
    const int spline_npoints_ind = source_species.findAttribute("spline_npoints");
    const int lmax_ind           = source_species.findAttribute("lmax");

    for (int center_idx = 0; center_idx < AtomicCentersInfo.Ncenters; center_idx++)
    {
      const int my_GroupID = AtomicCentersInfo.GroupID[center_idx];
      if (inner_cutoff_ind >= 0)
        AtomicCentersInfo.inner_cutoff[center_idx] = source_species(inner_cutoff_ind, my_GroupID);
      if (cutoff_radius_ind >= 0)
        AtomicCentersInfo.cutoff[center_idx] = source_species(cutoff_radius_ind, my_GroupID);
      if (spline_radius_ind >= 0)
        AtomicCentersInfo.spline_radius[center_idx] = source_species(spline_radius_ind, my_GroupID);
      if (spline_npoints_ind >= 0)
        AtomicCentersInfo.spline_npoints[center_idx] = source_species(spline_npoints_ind, my_GroupID);
      if (lmax_ind >= 0)
        AtomicCentersInfo.lmax[center_idx] = source_species(lmax_ind, my_GroupID);
    }
  }
  ///////////////////////////
  // Read the twist angles //
  ///////////////////////////
  primcell_kpoints.resize(NumTwists);
  for (int ti = 0; ti < NumTwists; ti++)
  {
    std::ostringstream path;
    path << "/electrons/kpoint_" << ti << "/reduced_k";
    TinyVector<double, OHMMS_DIM> primcell_kpoints_DP;
    H5File.read(primcell_kpoints_DP, path.str());
    primcell_kpoints[ti] = primcell_kpoints_DP;
  }
  if (qmc_common.use_density)
  {
    //////////////////////////////////////////////////////////
    // Only if it is bulk: If the density has not been set in TargetPtcl, and   //
    // the density is available, read it in and save it     //
    // in TargetPtcl.                                       //
    //////////////////////////////////////////////////////////
    if (TargetPtcl.getLattice().SuperCellEnum == SUPERCELL_BULK)
    {
      // FIXME:  add support for more than one spin density
      if (TargetPtcl.Density_G.empty())
      {
        Array<double, OHMMS_DIM> Density_r_DP;
        TinyVector<int, 3> mesh;
        H5File.read(TargetPtcl.DensityReducedGvecs, "/electrons/density/gvectors");
        int numG = TargetPtcl.DensityReducedGvecs.size();
// Convert primitive G-vectors to supercell G-vectors
// Also, flip sign since ESHDF format uses opposite sign convention
#pragma omp parallel for
        for (int iG = 0; iG < numG; iG++)
          TargetPtcl.DensityReducedGvecs[iG] = -1 * dot(TileMatrix, TargetPtcl.DensityReducedGvecs[iG]);
        app_log() << "  Read " << numG << " density G-vectors.\n";
        for (int ispin = 0; ispin < NumSpins; ispin++)
        {
          std::ostringstream density_r_path, density_g_path;
          density_r_path << "/electrons/density/spin_" << ispin << "/density_r";
          density_g_path << "/electrons/density/spin_" << ispin << "/density_g";
          H5File.readEntry(Density_r_DP, density_r_path.str());
          TargetPtcl.Density_r = Density_r_DP;
          if (TargetPtcl.DensityReducedGvecs.size())
          {
            app_log() << "  EinsplineSetBuilder found density in the HDF5 file.\n";
            std::vector<ComplexType> density_G;
            std::vector<std::complex<double>> Density_G_DP;
            H5File.read(Density_G_DP, density_g_path.str());
            density_G.assign(Density_G_DP.begin(), Density_G_DP.end());
            if (!density_G.size())
            {
              app_error() << "  Density reduced G-vectors defined, but not the"
                          << " density.\n";
              abort();
            }
            else
            {
              if (ispin == 0)
                TargetPtcl.Density_G = density_G;
              else
                for (int iG = 0; iG < density_G.size(); iG++)
                  TargetPtcl.Density_G[iG] += density_G[iG];
            }
          }
        }
      }
      //////////////////////////////////////////////////////////
      // If the density has not been set in TargetPtcl, and   //
      // the density is available, read it in and save it     //
      // in TargetPtcl.                                       //
      //////////////////////////////////////////////////////////
      // FIXME:  add support for more than one spin potential
      if (!TargetPtcl.VHXC_r[0].size())
      {
        TinyVector<int, 3> mesh;
        H5File.readEntry(TargetPtcl.VHXCReducedGvecs, "/electrons/VHXC/gvectors");
        int numG = TargetPtcl.VHXCReducedGvecs.size();
// Convert primitive G-vectors to supercell G-vectors
// Also, flip sign since ESHDF format uses opposite sign convention
#pragma omp parallel for
        for (int iG = 0; iG < numG; iG++)
          TargetPtcl.VHXCReducedGvecs[iG] = -1 * dot(TileMatrix, TargetPtcl.VHXCReducedGvecs[iG]);
        app_log() << "  Read " << numG << " VHXC G-vectors.\n";
        for (int ispin = 0; ispin < NumSpins; ispin++)
        {
          Array<double, OHMMS_DIM> VHXC_r_DP;
          std::ostringstream VHXC_r_path, VHXC_g_path;
          VHXC_r_path << "/electrons/VHXC/spin_" << ispin << "/VHXC_r";
          VHXC_g_path << "/electrons/VHXC/spin_" << ispin << "/VHXC_g";
          H5File.readEntry(VHXC_r_DP, VHXC_r_path.str());
          TargetPtcl.VHXC_r[ispin] = VHXC_r_DP;
          if (TargetPtcl.VHXCReducedGvecs.size())
          {
            app_log() << "  EinsplineSetBuilder found VHXC in the HDF5 file.\n";
            std::vector<std::complex<double>> VHXC_G_DP;
            std::vector<ComplexType> VHXC_G;
            H5File.read(VHXC_G_DP, VHXC_g_path.str());
            VHXC_G.assign(VHXC_G_DP.begin(), VHXC_G_DP.end());
            if (!VHXC_G.size())
            {
              app_error() << "  VHXC reduced G-vectors defined, but not the"
                          << " VHXC.\n";
              abort();
            }
            else
              TargetPtcl.VHXC_G[ispin] = VHXC_G;
          }
        }
      }
    }
  }
  else
  {
    app_log() << "   Skip initialization of the density" << std::endl;
  }
  return true;
}

bool EinsplineSetBuilder::ReadGvectors_ESHDF()
{
  bool root = myComm->rank() == 0;
  //this is always ugly
  MeshSize    = 0;
  int hasPsig = 1;
  if (root)
  {
    H5File.readEntry(MeshSize, "/electrons/psi_r_mesh");
    H5File.readEntry(MeshSize, "/electrons/mesh");
  }
  myComm->bcast(MeshSize);
  hasPsig = (MeshSize[0] == 0);
  if (hasPsig)
  {
    int nallowed  = 257;
    int allowed[] = {72,    75,    80,    81,    90,    96,    100,   108,   120,   125,   128,   135,   144,   150,
                     160,   162,   180,   192,   200,   216,   225,   240,   243,   250,   256,   270,   288,   300,
                     320,   324,   360,   375,   384,   400,   405,   432,   450,   480,   486,   500,   512,   540,
                     576,   600,   625,   640,   648,   675,   720,   729,   750,   768,   800,   810,   864,   900,
                     960,   972,   1000,  1024,  1080,  1125,  1152,  1200,  1215,  1250,  1280,  1296,  1350,  1440,
                     1458,  1500,  1536,  1600,  1620,  1728,  1800,  1875,  1920,  1944,  2000,  2025,  2048,  2160,
                     2187,  2250,  2304,  2400,  2430,  2500,  2560,  2592,  2700,  2880,  2916,  3000,  3072,  3125,
                     3200,  3240,  3375,  3456,  3600,  3645,  3750,  3840,  3888,  4000,  4050,  4096,  4320,  4374,
                     4500,  4608,  4800,  4860,  5000,  5120,  5184,  5400,  5625,  5760,  5832,  6000,  6075,  6144,
                     6250,  6400,  6480,  6561,  6750,  6912,  7200,  7290,  7500,  7680,  7776,  8000,  8100,  8192,
                     8640,  8748,  9000,  9216,  9375,  9600,  9720,  10000, 10125, 10240, 10368, 10800, 10935, 11250,
                     11520, 11664, 12000, 12150, 12288, 12500, 12800, 12960, 13122, 13500, 13824, 14400, 14580, 15000,
                     15360, 15552, 15625, 16000, 16200, 16384, 16875, 17280, 17496, 18000, 18225, 18432, 18750, 19200,
                     19440, 19683, 20000, 20250, 20480, 20736, 21600, 21870, 22500, 23040, 23328, 24000, 24300, 24576,
                     25000, 25600, 25920, 26244, 27000, 27648, 28125, 28800, 29160, 30000, 30375, 30720, 31104, 31250,
                     32000, 32400, 32768, 32805, 33750, 34560, 34992, 36000, 36450, 36864, 37500, 38400, 38880, 39366,
                     40000, 40500, 40960, 41472, 43200, 43740, 45000, 46080, 46656, 46875, 48000, 48600, 49152, 50000,
                     50625, 51200, 51840, 52488, 54000, 54675, 55296, 56250, 57600, 58320, 59049, 60000, 60750, 61440,
                     62208, 62500, 64000, 64800, 65536};
    MaxNumGvecs   = 0;
    //    std::set<TinyVector<int,3> > Gset;
    // Read k-points for all G-vectors and take the union
    TinyVector<int, 3> maxIndex(0, 0, 0);
    Gvecs.resize(NumTwists);
    {
      int numg = 0;
      if (root)
      {
        std::ostringstream Gpath;
        Gpath << "/electrons/kpoint_0/gvectors";
        H5File.read(Gvecs[0], Gpath.str());
        numg = Gvecs[0].size();
      }
      myComm->bcast(numg);
      if (!root)
        Gvecs[0].resize(numg);
      myComm->bcast(Gvecs[0]);
      MaxNumGvecs = Gvecs[0].size();
      for (int ig = 0; ig < Gvecs[0].size(); ig++)
      {
        maxIndex[0] = std::max(maxIndex[0], std::abs(Gvecs[0][ig][0]));
        maxIndex[1] = std::max(maxIndex[1], std::abs(Gvecs[0][ig][1]));
        maxIndex[2] = std::max(maxIndex[2], std::abs(Gvecs[0][ig][2]));
      }
      // for (int ig=0; ig<Gvecs.size(); ig++)
      // 	if (Gset.find(Gvecs[ig]) == Gset.end())
      // 	  Gset.insert(Gvecs[ig]);
    } //done with kpoint_0
    MeshSize[0] = (int)std::ceil(4.0 * MeshFactor * maxIndex[0]);
    MeshSize[1] = (int)std::ceil(4.0 * MeshFactor * maxIndex[1]);
    MeshSize[2] = (int)std::ceil(4.0 * MeshFactor * maxIndex[2]);
    //only use 2^a 3^b 5^c where a>=2  up to 65536
    int* ix     = std::lower_bound(allowed, allowed + nallowed, MeshSize[0]);
    int* iy     = std::lower_bound(allowed, allowed + nallowed, MeshSize[1]);
    int* iz     = std::lower_bound(allowed, allowed + nallowed, MeshSize[2]);
    MeshSize[0] = (MeshSize[0] > 128) ? *ix : (MeshSize[0] + MeshSize[0] % 2);
    MeshSize[1] = (MeshSize[1] > 128) ? *iy : (MeshSize[1] + MeshSize[1] % 2);
    MeshSize[2] = (MeshSize[2] > 128) ? *iz : (MeshSize[2] + MeshSize[2] % 2);
    if (Version[0] < 2)
    {
      //get the map for each twist, but use the MeshSize from kpoint_0
      app_log() << "  ESHDF::Version " << Version << std::endl;
      app_log() << "  Assumes distinct Gvecs set for different twists. Regenerate orbital files using updated QE."
                << std::endl;
      for (int k = 0; k < DistinctTwists.size(); ++k)
      {
        int ik = DistinctTwists[k];
        if (ik == 0)
          continue; //already done
        int numg = 0;
        if (root)
        {
          std::ostringstream Gpath;
          Gpath << "/electrons/kpoint_" << ik << "/gvectors";
          H5File.read(Gvecs[ik], Gpath.str());
          numg = Gvecs[ik].size();
        }
        myComm->bcast(numg);
        if (numg == 0)
        {
          //copy kpoint_0, default
          Gvecs[ik] = Gvecs[0];
        }
        else
        {
          if (numg != MaxNumGvecs)
          {
            std::ostringstream o;
            o << "Twist " << ik << ": The number of Gvecs is different from kpoint_0."
              << " This is not supported anymore. Rerun pw2qmcpack.x or equivalent";
            APP_ABORT(o.str());
          }
          if (!root)
            Gvecs[ik].resize(numg);
          myComm->bcast(Gvecs[ik]);
        }
      }
    }
  }
  app_log() << "B-spline mesh factor is " << MeshFactor << std::endl;
  app_log() << "B-spline mesh size is (" << MeshSize[0] << ", " << MeshSize[1] << ", " << MeshSize[2] << ")\n";
  app_log() << "Maxmimum number of Gvecs " << MaxNumGvecs << std::endl;
  app_log().flush();
  return hasPsig;
}

void EinsplineSetBuilder::OccupyBands_ESHDF(int spin, int sortBands, int numOrbs)
{
  if (myComm->rank() != 0)
    return;

  std::vector<BandInfo>& SortBands(*FullBands[spin]);
  SortBands.clear(); //??? can exit if SortBands is already made?
  int maxOrbs(0);
  for (int ti = 0; ti < DistinctTwists.size(); ti++)
  {
    int tindex = DistinctTwists[ti];
    // First, read valence states
    std::ostringstream ePath;
    ePath << "/electrons/kpoint_" << tindex << "/spin_" << spin << "/eigenvalues";
    std::vector<double> eigvals;
    H5File.read(eigvals, ePath.str());
    for (int bi = 0; bi < NumBands; bi++)
    {
      BandInfo band;
      band.TwistIndex    = tindex;
      band.BandIndex     = bi;
      band.MakeTwoCopies = MakeTwoCopies[ti];
      band.Energy        = eigvals[bi];
      if (band.Energy > -1.0e100)
        SortBands.push_back(band);
      if (MakeTwoCopies[ti])
        maxOrbs += 2;
      else
        maxOrbs++;
    }
  }

  app_log() << SortBands.size() << " complex-valued orbitals supplied by h5 can be expanded up to " << maxOrbs
            << " SPOs." << std::endl;
  if (maxOrbs < numOrbs)
    myComm->barrier_and_abort("EinsplineSetBuilder::OccupyBands_ESHDF user input requests "
                              "more orbitals than what the h5 file supplies.");

  // Now sort the bands by energy
  if (sortBands == 2)
  {
    app_log() << "Sorting the bands by index now:\n";
    sort(SortBands.begin(), SortBands.end(), sortByIndex);
  }
  else if (sortBands == 1)
  {
    app_log() << "Sorting the bands now:\n";
    sort(SortBands.begin(), SortBands.end());
  }

  std::vector<int> gsOcc(maxOrbs);
  int N_gs_orbs = numOrbs;
  int nocced(0);
  for (int ti = 0; ti < SortBands.size(); ti++)
  {
    if (nocced < N_gs_orbs)
    {
      if (SortBands[ti].MakeTwoCopies && (N_gs_orbs - nocced > 1))
      {
        nocced += 2;
        gsOcc[ti] = 2;
      }
      else if ((SortBands[ti].MakeTwoCopies && (N_gs_orbs - nocced == 1)) || !SortBands[ti].MakeTwoCopies)
      {
        nocced += 1;
        gsOcc[ti] = 1;
      }
    }
  }
  if (occ_format == "energy")
  {
    app_log() << "  Occupying bands based on energy in mode " << (Occ.size() > 0 ? "\"excited\"" : "\"ground\"")
              << std::endl;
    // To get the occupations right.
    std::vector<int> Removed(0, 0);
    std::vector<int> Added(0, 0);
    for (int ien = 0; ien < Occ.size(); ien++)
    {
      if (Occ[ien] < 0)
        Removed.push_back(-Occ[ien]);
      else if (Occ[ien] > 0)
        Added.push_back(Occ[ien]);
    }
    if (Added.size() - Removed.size() != 0)
    {
      app_log() << "need to add and remove same number of orbitals. " << Added.size() << " " << Removed.size()
                << std::endl;
      APP_ABORT("ChangedOccupations");
    }
    std::vector<int> DiffOcc(maxOrbs, 0);
    //Probably a cleaner way to do this.
    for (int i = 0; i < Removed.size(); i++)
      DiffOcc[Removed[i] - 1] -= 1;
    for (int i = 0; i < Added.size(); i++)
      DiffOcc[Added[i] - 1] += 1;
    std::vector<int> SumOrb(SortBands.size(), 0);
    int doi(0);
    for (int i = 0; i < SumOrb.size(); i++)
    {
      if (SortBands[i].MakeTwoCopies)
      {
        SumOrb[i] = gsOcc[i] + DiffOcc[doi++];
        SumOrb[i] += DiffOcc[doi++];
      }
      else
        SumOrb[i] = gsOcc[i] + DiffOcc[doi++];
    }
    std::vector<BandInfo> ReOrderedBands;
    std::vector<BandInfo> RejectedBands;
    for (int i = 0; i < SumOrb.size(); i++)
    {
      if (SumOrb[i] == 2)
      {
        SortBands[i].MakeTwoCopies = true;
        ReOrderedBands.push_back(SortBands[i]);
      }
      else if (SumOrb[i] == 1)
      {
        SortBands[i].MakeTwoCopies = false;
        ReOrderedBands.push_back(SortBands[i]);
      }
      else if (SumOrb[i] == 0)
      {
        SortBands[i].MakeTwoCopies = false;
        RejectedBands.push_back(SortBands[i]);
      }
      else
      {
        app_log() << " Trying to add the same orbital (" << i << ") less than zero or more than 2 times." << std::endl;
        APP_ABORT("Sorting Excitation");
      }
    }
    ReOrderedBands.insert(ReOrderedBands.end(), RejectedBands.begin(), RejectedBands.end());
    SortBands = ReOrderedBands;
  }
  else if (occ_format == "band")
  {
    app_log() << "  Occupying bands based on (ti,bi) data." << std::endl;
    if (Occ.size() != particle_hole_pairs * 4)
    {
      app_log() << " Need Occ = pairs*4. Occ is (ti,bi) of removed, then added." << std::endl;
      app_log() << Occ.size() << " " << particle_hole_pairs << std::endl;
      APP_ABORT("ChangedOccupations");
    }
    int cnt(0);
    for (int ien = 0; ien < SortBands.size(); ien++)
    {
      if ((Occ[cnt] == SortBands[ien].TwistIndex) && (Occ[cnt + 1] == SortBands[ien].BandIndex))
      {
        if (cnt < particle_hole_pairs * 2)
        {
          gsOcc[ien] -= 1;
          cnt += 2;
          app_log() << "removing orbital " << ien << std::endl;
        }
        else
        {
          gsOcc[ien] += 1;
          app_log() << "adding orbital " << ien << std::endl;
          cnt += 2;
        }
      }
    }
    std::vector<BandInfo> ReOrderedBands;
    std::vector<BandInfo> RejectedBands;
    for (int i = 0; i < SortBands.size(); i++)
    {
      if (gsOcc[i] == 2)
      {
        SortBands[i].MakeTwoCopies = true;
        ReOrderedBands.push_back(SortBands[i]);
      }
      else if (gsOcc[i] == 1)
      {
        SortBands[i].MakeTwoCopies = false;
        ReOrderedBands.push_back(SortBands[i]);
      }
      else if (gsOcc[i] == 0)
      {
        SortBands[i].MakeTwoCopies = false;
        RejectedBands.push_back(SortBands[i]);
      }
      else
      {
        app_log() << " Trying to add the same orbital (" << i << ") less than zero or more than 2 times." << std::endl;
        APP_ABORT("Sorting Excitation");
      }
    }
    ReOrderedBands.insert(ReOrderedBands.end(), RejectedBands.begin(), RejectedBands.end());
    SortBands = ReOrderedBands;
  }
  //for(int sw=0;sw<Removed.size();sw++){
  //  app_log()<<" Swapping two orbitals "<<Removed[sw]<<" and "<<Added[sw]<< std::endl;
  //  BandInfo tempband(SortBands[Removed[sw]-1]);
  //  SortBands[Removed[sw]-1] = SortBands[Added[sw]-1];
  //  SortBands[Added[sw]-1] = tempband;
  //}
  int orbIndex        = 0;
  int numOrbs_counter = 0;
  while (numOrbs_counter < numOrbs)
  {
    if (SortBands[orbIndex].MakeTwoCopies)
      numOrbs_counter += 2;
    else
      numOrbs_counter++;
    orbIndex++;
  }
  NumDistinctOrbitals = orbIndex;
  app_log() << "We will read " << NumDistinctOrbitals << " distinct complex-valued orbitals from h5.\n";
}

} // namespace qmcplusplus
