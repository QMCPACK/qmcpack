//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Paul R. C. Kent, kentpr@ornl.gov, Oak Ridge National Laboratory
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Ying Wai Li, yingwaili@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "EinsplineSetBuilder.h"
#include <PlatformSelector.hpp>
#include "CPU/e2iphi.h"
#include "CPU/SIMD/vmath.hpp"
#include "OhmmsData/AttributeSet.h"
#include "Message/CommOperators.h"
#include "Utilities/Timer.h"
#include "ParticleBase/RandomSeqGenerator.h"
#include "Particle/DistanceTable.h"
#include "einspline_helper.hpp"
#include "BsplineReader.h"
#include "BsplineSet.h"
#include "createBsplineReader.h"

#include <array>
#include <string_view>

namespace qmcplusplus
{
void EinsplineSetBuilder::set_metadata(int numOrbs,
                                       int twist_num_inp,
                                       const TinyVector<double, OHMMS_DIM>& twist_inp,
                                       bool skipChecks)
{
  // 1. set a lot of internal parameters in the EinsplineSetBuilder class
  //  e.g. TileMatrix, use_real_splines_, DistinctTwists, MakeTwoCopies.
  // 2. this is also where metadata for the orbitals are read from the wavefunction hdf5 file
  //  and broadcast to MPI groups. Variables broadcasted are listed in
  //  EinsplineSetBuilderCommon.cpp EinsplineSetBuilder::BroadcastOrbitalInfo()
  //

  Timer orb_info_timer;
  // The tiling can be set by a simple vector, (e.g. 2x2x2), or by a
  // full 3x3 matrix of integers.  If the tilematrix was not set in
  // the input file...
  bool matrixNotSet = true;
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      matrixNotSet = matrixNotSet && (TileMatrix(i, j) == 0);
  // then set the matrix to identity.
  if (matrixNotSet)
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        TileMatrix(i, j) = (i == j) ? 1 : 0;
  if (myComm->rank() == 0)
  {
    std::array<char, 1000> buff;
    int length =
        std::snprintf(buff.data(), buff.size(), "  TileMatrix = \n [ %2d %2d %2d\n   %2d %2d %2d\n   %2d %2d %2d ]\n",
                      TileMatrix(0, 0), TileMatrix(0, 1), TileMatrix(0, 2), TileMatrix(1, 0), TileMatrix(1, 1),
                      TileMatrix(1, 2), TileMatrix(2, 0), TileMatrix(2, 1), TileMatrix(2, 2));
    if (length < 0)
      throw std::runtime_error("Error converting TileMatrix to a string");
    app_log() << std::string_view(buff.data(), length);
  }
  if (numOrbs == 0)
    myComm->barrier_and_abort(
        "EinsplineSetBuilder::createSPOSet You must specify the number of orbitals in the input file.");
  else
    app_log() << "  Reading " << numOrbs << " orbitals from HDF5 file.\n";

  /////////////////////////////////////////////////////////////////
  // Read the basic orbital information, without reading all the //
  // orbitals themselves.                                        //
  /////////////////////////////////////////////////////////////////
  orb_info_timer.restart();
  if (myComm->rank() == 0)
    if (!ReadOrbitalInfo(skipChecks))
      throw std::runtime_error("EinsplineSetBuilder::set_metadata Error reading orbital info from HDF5 file.");
  app_log() << "TIMER  EinsplineSetBuilder::ReadOrbitalInfo " << orb_info_timer.elapsed() << std::endl;
  myComm->barrier();

  orb_info_timer.restart();
  BroadcastOrbitalInfo();
  app_log() << "TIMER  EinsplineSetBuilder::BroadcastOrbitalInfo " << orb_info_timer.elapsed() << std::endl;
  app_log().flush();

  // setup primitive cell and supercell
  PrimCell.set(Lattice);
  SuperCell.set(SuperLattice);
  GGt = dot(transpose(PrimCell.G), PrimCell.G);

  // Now, analyze the k-point mesh to figure out the what k-points  are needed
  AnalyzeTwists2(twist_num_inp, twist_inp);
}

std::unique_ptr<SPOSet> EinsplineSetBuilder::createSPOSetFromXML(xmlNodePtr cur)
{
  //use 2 bohr as the default when truncated orbitals are used based on the extend of the ions
  int numOrbs = 0;
  int sortBands(1);
  int spinSet       = 0;
  bool skipChecks   = false;
  int twist_num_inp = TWISTNUM_NO_INPUT;
  TinyVector<double, OHMMS_DIM> twist_inp(TWIST_NO_INPUT);

  std::string sourceName;
  std::string spo_prec("double");
  std::string truncate("no");
  std::string hybrid_rep("no");
  std::string skip_checks("no");
  std::string use_einspline_set_extended(
      "no"); // use old spline library for high-order derivatives, e.g. needed for backflow optimization
  std::string useGPU;
  std::string GPUsharing = "no";

  ScopedTimer spo_timer_scope(createGlobalTimer("einspline::CreateSPOSetFromXML", timer_level_medium));

  {
    TinyVector<int, OHMMS_DIM> TileFactor_do_not_use;
    OhmmsAttributeSet a;
    a.add(H5FileName, "href");
    a.add(TileFactor_do_not_use, "tile", {}, TagStatus::DELETED);
    a.add(sortBands, "sort");
    a.add(TileMatrix, "tilematrix");
    a.add(twist_num_inp, "twistnum");
    a.add(twist_inp, "twist");
    a.add(sourceName, "source");
    a.add(MeshFactor, "meshfactor");
    a.add(hybrid_rep, "hybridrep");
    a.add(useGPU, "gpu", CPUOMPTargetSelector::candidate_values);
    a.add(GPUsharing, "gpusharing"); // split spline across GPUs visible per rank
    a.add(spo_prec, "precision");
    a.add(truncate, "truncate");
    a.add(myName, "tag");
    a.add(skip_checks, "skip_checks");

    a.put(XMLRoot);
    a.add(numOrbs, "size");
    a.add(numOrbs, "norbs");
    a.add(spinSet, "spindataset");
    a.add(spinSet, "group");
    a.put(cur);

    if (myName.empty())
      myName = "einspline";
  }

  if (skip_checks == "yes")
    skipChecks = true;

  auto pit(ParticleSets.find(sourceName));
  if (pit == ParticleSets.end())
    myComm->barrier_and_abort("Einspline needs the source particleset");
  else
    SourcePtcl = pit->second.get();

  ///////////////////////////////////////////////
  // Read occupation information from XML file //
  ///////////////////////////////////////////////
  const std::vector<int> last_occ(Occ);
  Occ.resize(0, 0); // correspond to ground
  bool NewOcc(false);

  {
    OhmmsAttributeSet oAttrib;
    oAttrib.add(spinSet, "spindataset");
    oAttrib.put(cur);
  }

  xmlNodePtr spo_cur = cur;
  cur                = cur->children;
  while (cur != NULL)
  {
    std::string cname((const char*)(cur->name));
    if (cname == "occupation")
    {
      std::string occ_mode("ground");
      occ_format          = "energy";
      particle_hole_pairs = 0;
      OhmmsAttributeSet oAttrib;
      oAttrib.add(occ_mode, "mode");
      oAttrib.add(spinSet, "spindataset");
      oAttrib.add(occ_format, "format");
      oAttrib.add(particle_hole_pairs, "pairs");
      oAttrib.put(cur);
      if (occ_mode == "excited")
        putContent(Occ, cur);
      else if (occ_mode != "ground")
        myComm->barrier_and_abort("EinsplineSetBuilder::createSPOSet Only ground state occupation "
                                  "currently supported in EinsplineSetBuilder.");
    }
    cur = cur->next;
  }
  if (Occ != last_occ)
  {
    NewOcc = true;
  }
  else
    NewOcc = false;
#if defined(MIXED_PRECISION)
  app_log() << "\t  MIXED_PRECISION=1 Overwriting the einspline storage to single precision.\n";
  spo_prec = "single"; //overwrite
#endif
  H5OrbSet aset(H5FileName, spinSet, numOrbs);
  const auto iter = SPOSetMap.find(aset);
  if ((iter != SPOSetMap.end()) && (!NewOcc))
    app_warning() << "!!!!!!! Identical SPOSets are detected by EinsplineSetBuilder! "
                     "Implicit sharing one SPOSet for spin-up and spin-down electrons has been removed. "
                     "Each determinant creates its own SPOSet with dedicated memory for spline coefficients. "
                     "To avoid increasing the memory footprint of spline coefficients, "
                     "create a single SPOset outside the determinantset using 'sposet_collection' "
                     "and reference it by name on the determinant line."
                  << std::endl;

  if (FullBands[spinSet] == 0)
    FullBands[spinSet] = std::make_unique<std::vector<BandInfo>>();

  // Ensure the first SPO set must be spinSet==0
  // to correctly initialize key data of EinsplineSetBuilder
  if (SPOSetMap.size() == 0 && spinSet != 0)
    myComm->barrier_and_abort("The first SPO set must have spindataset=\"0\"");

  // set the internal parameters
  if (spinSet == 0)
    set_metadata(numOrbs, twist_num_inp, twist_inp, skipChecks);

  //////////////////////////////////
  // Create the OrbitalSet object
  //////////////////////////////////
  Timer mytimer;
  mytimer.restart();
  OccupyBands(spinSet, sortBands, numOrbs, skipChecks);
  if (spinSet == 0)
    TileIons();

  bool use_single = (spo_prec == "single" || spo_prec == "float");

  // safeguard for a removed feature
  if (truncate == "yes")
    myComm->barrier_and_abort(
        "The 'truncate' feature of spline SPO has been removed. Please use hybrid orbital representation.");

#if !defined(QMC_COMPLEX)
  if (use_real_splines_)
  {
    //if(TargetPtcl.Lattice.SuperCellEnum != SUPERCELL_BULK && truncate=="yes")
    if (MixedSplineReader == 0)
      MixedSplineReader = createBsplineReal(this, use_single, hybrid_rep == "yes", useGPU);
  }
  else
#endif
  {
    if (MixedSplineReader == 0)
      MixedSplineReader = createBsplineComplex(this, use_single, hybrid_rep == "yes", useGPU);
  }

  MixedSplineReader->setCommon(XMLRoot);
  // temporary disable the following function call, Ye Luo
  // RotateBands_ESHDF(spinSet, dynamic_cast<EinsplineSetExtended<std::complex<double> >*>(OrbitalSet));
  bcastSortBands(spinSet, NumDistinctOrbitals, myComm->rank() == 0);
  auto OrbitalSet = MixedSplineReader->create_spline_set(spinSet, spo_cur);
  if (!OrbitalSet)
    myComm->barrier_and_abort("Failed to create SPOSet*");
  app_log() << "Time spent in creating B-spline SPOs " << mytimer.elapsed() << " sec" << std::endl;
  OrbitalSet->finalizeConstruction();
  SPOSetMap[aset] = OrbitalSet.get();
  return OrbitalSet;
}

std::unique_ptr<SPOSet> EinsplineSetBuilder::createSPOSet(xmlNodePtr cur, SPOSetInputInfo& input_info)
{
  if (MixedSplineReader == 0)
    myComm->barrier_and_abort("EinsplineSetExtended<T> cannot create a SPOSet");

  std::string aname;
  int spinSet(0);
  OhmmsAttributeSet a;
  a.add(spinSet, "spindataset");
  a.add(spinSet, "group");
  a.put(cur);

  //allow only non-overlapping index sets and use the max index as the identifier
  int norb = input_info.max_index();
  H5OrbSet aset(H5FileName, spinSet, norb);

  auto bspline_zd = MixedSplineReader->create_spline_set(spinSet, cur, input_info);
  if (bspline_zd)
    SPOSetMap[aset] = bspline_zd.get();
  return bspline_zd;
}

} // namespace qmcplusplus
