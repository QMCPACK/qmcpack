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


#include "QMCWaveFunctions/EinsplineSetBuilder.h"
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
  int have_dpsi         = false;
  int NumAtomicOrbitals = 0;
  NumCoreStates = NumMuffinTins = NumTwists = NumSpins = NumBands = NumAtomicOrbitals = 0;
  NumElectrons = TargetPtcl.getTotalNum();
  H5File.read(NumBands, "/electrons/kpoint_0/spin_0/number_of_states");
  H5File.readEntry(NumCoreStates, "/electrons/kpoint_0/spin_0/number_of_core_states");
  H5File.readEntry(NumSpins, "/electrons/number_of_spins");
  H5File.read(NumTwists, "/electrons/number_of_kpoints");
  H5File.readEntry(NumMuffinTins, "/muffin_tins/number_of_tins");
  H5File.readEntry(have_dpsi, "/electrons/have_dpsi");
  H5File.readEntry(NumAtomicOrbitals, "/electrons/number_of_atomic_orbitals");
  HaveOrbDerivs = have_dpsi;
  app_log() << "bands=" << NumBands << ", elecs=" << NumElectrons << ", spins=" << NumSpins << ", twists=" << NumTwists
            << ", muffin tins=" << NumMuffinTins << ", core states=" << NumCoreStates << std::endl;
  app_log() << "atomic orbital=" << NumAtomicOrbitals << std::endl;
  if (TileFactor[0] != 1 || TileFactor[1] != 1 || TileFactor[2] != 1)
    app_log() << "  Using a " << TileFactor[0] << "x" << TileFactor[1] << "x" << TileFactor[2] << " tiling factor.\n";
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
  /////////////////////////////////////
  // Read atomic orbital information //
  /////////////////////////////////////
  AtomicOrbitals.resize(NumAtomicOrbitals);
  for (int iat = 0; iat < NumAtomicOrbitals; iat++)
  {
    AtomicOrbital<std::complex<double>>& orb = AtomicOrbitals[iat];
    int lmax, polynomial_order, spline_points;
    RealType cutoff_radius, polynomial_radius, spline_radius;
    PosType position;
    double cutoff_radius_DP, polynomial_radius_DP, spline_radius_DP;
    TinyVector<double, OHMMS_DIM> position_DP;
    std::ostringstream groupstream;
    groupstream << "/electrons/atomic_orbital_" << iat << "/";
    std::string groupname = groupstream.str();
    H5File.read(lmax, groupname + "lmax");
    H5File.read(polynomial_order, groupname + "polynomial_order");
    H5File.read(spline_points, groupname + "spline_points");
    H5File.read(cutoff_radius_DP, groupname + "cutoff_radius");
    H5File.read(polynomial_radius_DP, groupname + "polynomial_radius");
    H5File.read(spline_radius_DP, groupname + "spline_radius");
    H5File.read(position_DP, groupname + "position");
    cutoff_radius     = cutoff_radius_DP;
    polynomial_radius = polynomial_radius_DP;
    spline_radius     = spline_radius_DP;
    position          = position_DP;
    orb.set_pos(position);
    orb.set_lmax(lmax);
    orb.set_cutoff(cutoff_radius);
    orb.set_spline(spline_radius, spline_points);
    orb.set_polynomial(polynomial_radius, polynomial_order);
  }
  ///////////////////////////
  // Read the twist angles //
  ///////////////////////////
  TwistAngles.resize(NumTwists);
  TwistSymmetry.resize(NumTwists);
  TwistWeight.resize(NumTwists);
  for (int ti = 0; ti < NumTwists; ti++)
  {
    std::ostringstream path;
    path << "/electrons/kpoint_" << ti << "/reduced_k";
    TinyVector<double, OHMMS_DIM> TwistAngles_DP;
    H5File.read(TwistAngles_DP, path.str());
    TwistAngles[ti] = TwistAngles_DP;
    if ((Version[0] >= 2) and (Version[1] >= 1))
    {
      std::ostringstream sym_path;
      sym_path << "/electrons/kpoint_" << ti << "/symgroup";
      H5File.readEntry(TwistSymmetry[ti], sym_path.str());
      std::ostringstream nsym_path;
      nsym_path << "/electrons/kpoint_" << ti << "/numsym";
      H5File.readEntry(TwistWeight[ti], nsym_path.str());
    }
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
      band.IsCoreState   = false;
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
    // Now, read core states
    for (int cs = 0; cs < NumCoreStates; cs++)
    {
      BandInfo band;
      band.IsCoreState   = true;
      band.TwistIndex    = tindex;
      band.BandIndex     = cs;
      band.MakeTwoCopies = MakeTwoCopies[ti];
      H5File.read(band.Energy, CoreStatePath(ti, cs) + "eigenvalue");
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
    app_log() << "  Occupying bands based on energy in mode " << (Occ.size() > 0? "\"excited\"" : "\"ground\"") << std::endl;
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
  NumValenceOrbs      = 0;
  NumCoreOrbs         = 0;
  while (numOrbs_counter < numOrbs)
  {
    if (SortBands[orbIndex].MakeTwoCopies)
      numOrbs_counter += 2;
    else
      numOrbs_counter++;
    if (SortBands[orbIndex].IsCoreState)
      NumCoreOrbs++;
    else
      NumValenceOrbs++;
    orbIndex++;
  }
  NumDistinctOrbitals = orbIndex;
  app_log() << "We will read " << NumDistinctOrbitals << " distinct complex-valued orbitals from h5.\n";
  app_log() << "There are " << NumCoreOrbs << " core states and " << NumValenceOrbs << " valence states.\n";
}

} // namespace qmcplusplus
