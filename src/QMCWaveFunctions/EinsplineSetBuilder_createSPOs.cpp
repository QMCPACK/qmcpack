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


#include "QMCWaveFunctions/EinsplineSetBuilder.h"
#include <PlatformSelector.hpp>
#include "CPU/e2iphi.h"
#include "CPU/SIMD/vmath.hpp"
#include "OhmmsData/AttributeSet.h"
#include "Message/CommOperators.h"
#include "Utilities/Timer.h"
#include "ParticleBase/RandomSeqGenerator.h"
#include "Particle/DistanceTable.h"
#include <fftw3.h>
#include "Utilities/ProgressReportEngine.h"
#include "QMCWaveFunctions/einspline_helper.hpp"
#if !defined(MIXED_PRECISION)
#include "QMCWaveFunctions/EinsplineSet.h"
#endif
#include "QMCWaveFunctions/BsplineFactory/BsplineReaderBase.h"
#include "QMCWaveFunctions/BsplineFactory/BsplineSet.h"
#include "QMCWaveFunctions/BsplineFactory/createBsplineReader.h"

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
  // then set the matrix to what may have been specified in the
  // tiling vector
  if (matrixNotSet)
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        TileMatrix(i, j) = (i == j) ? TileFactor[i] : 0;
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
  for (int iat = 0; iat < AtomicOrbitals.size(); iat++)
    AtomicOrbitals[iat].Lattice.set(Lattice);

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
  std::string spo_object_name;

  ScopedTimer spo_timer_scope(createGlobalTimer("einspline::CreateSPOSetFromXML", timer_level_medium));

  {
    OhmmsAttributeSet a;
    a.add(H5FileName, "href");
    a.add(TileFactor, "tile");
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
    a.add(use_einspline_set_extended, "use_old_spline");
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
    oAttrib.add(spo_object_name, "name");
    oAttrib.add(spo_object_name, "id");
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
  //if (use_complex_orb == "yes") use_real_splines_ = false; // override given user input

  // look for <backflow>, would be a lot easier with xpath, but I cannot get it to work
  bool has_backflow = false;

  xmlNodePtr wf  = XMLRoot->parent; // <wavefuntion>
  xmlNodePtr kid = wf->children;
  while (kid != NULL)
  {
    std::string tag((const char*)(kid->name));
    if (tag == "determinantset" || tag == "sposet_builder")
    {
      xmlNodePtr kid1 = kid->children;
      while (kid1 != NULL)
      {
        std::string tag1((const char*)(kid1->name));
        if (tag1 == "backflow")
        {
          has_backflow = true;
        }
        kid1 = kid1->next;
      }
    }
    kid = kid->next;
  }

  if (has_backflow && use_einspline_set_extended == "yes" && use_real_splines_)
    myComm->barrier_and_abort("backflow optimization is broken with use_real_splines_");

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
    {
      if (use_single)
        MixedSplineReader = createBsplineRealSingle(this, hybrid_rep == "yes", useGPU);
      else
        MixedSplineReader = createBsplineRealDouble(this, hybrid_rep == "yes", useGPU);
    }
  }
  else
#endif
  {
    if (MixedSplineReader == 0)
    {
      if (use_single)
        MixedSplineReader = createBsplineComplexSingle(this, hybrid_rep == "yes", useGPU);
      else
        MixedSplineReader = createBsplineComplexDouble(this, hybrid_rep == "yes", useGPU);
    }
  }

  MixedSplineReader->setCommon(XMLRoot);
  // temporary disable the following function call, Ye Luo
  // RotateBands_ESHDF(spinSet, dynamic_cast<EinsplineSetExtended<std::complex<double> >*>(OrbitalSet));
  HasCoreOrbs     = bcastSortBands(spinSet, NumDistinctOrbitals, myComm->rank() == 0);
  auto OrbitalSet = MixedSplineReader->create_spline_set(spinSet, spo_cur);
  if (!OrbitalSet)
    myComm->barrier_and_abort("Failed to create SPOSet*");
#if defined(MIXED_PRECISION)
  if (use_einspline_set_extended == "yes")
    myComm->barrier_and_abort("Option use_old_spline is not supported by the mixed precision build!");
#else
  if (use_einspline_set_extended == "yes")
  {
    std::unique_ptr<EinsplineSet> new_OrbitalSet;
    if (use_real_splines_)
    {
      auto temp_OrbitalSet                      = std::make_unique<EinsplineSetExtended<double>>(spo_object_name);
      temp_OrbitalSet->MultiSpline              = MixedSplineReader->export_MultiSplineDouble().release();
      temp_OrbitalSet->MultiSpline->num_splines = NumDistinctOrbitals;
      temp_OrbitalSet->resizeStorage(NumDistinctOrbitals, NumValenceOrbs);
      //set the flags for anti periodic boundary conditions
      temp_OrbitalSet->HalfG = dynamic_cast<BsplineSet&>(*OrbitalSet).getHalfG();
      new_OrbitalSet         = std::move(temp_OrbitalSet);
    }
    else
    {
      auto temp_OrbitalSet         = std::make_unique<EinsplineSetExtended<std::complex<double>>>(spo_object_name);
      temp_OrbitalSet->MultiSpline = MixedSplineReader->export_MultiSplineComplexDouble().release();
      temp_OrbitalSet->MultiSpline->num_splines = NumDistinctOrbitals;
      temp_OrbitalSet->resizeStorage(NumDistinctOrbitals, NumValenceOrbs);
      for (int iorb = 0, num = 0; iorb < NumDistinctOrbitals; iorb++)
      {
        int ti                               = (*FullBands[spinSet])[iorb].TwistIndex;
        temp_OrbitalSet->kPoints[iorb]       = PrimCell.k_cart(-TwistAngles[ti]);
        temp_OrbitalSet->MakeTwoCopies[iorb] = (num < (numOrbs - 1)) && (*FullBands[spinSet])[iorb].MakeTwoCopies;
        num += temp_OrbitalSet->MakeTwoCopies[iorb] ? 2 : 1;
      }
      new_OrbitalSet = std::move(temp_OrbitalSet);
    }
    //set the internal parameters
    setTiling(new_OrbitalSet.get(), numOrbs);
    OrbitalSet = std::move(new_OrbitalSet);
  }
#endif
  app_log() << "Time spent in creating B-spline SPOs " << mytimer.elapsed() << "sec" << std::endl;
#ifdef Ye_debug
#ifndef QMC_COMPLEX
  if (myComm->rank() == 0 && OrbitalSet->MuffinTins.size() > 0)
  {
    FILE* fout = fopen("TestMuffins.dat", "w");
    Vector<double> phi(numOrbs), lapl(numOrbs);
    Vector<PosType> grad(numOrbs);
    ParticleSet P;
    P.R.resize(6);
    for (int i = 0; i < P.R.size(); i++)
      P.R[i] = PosType(0.0, 0.0, 0.0);
    PosType N = 0.25 * PrimCell.a(0) + 0.25 * PrimCell.a(1) + 0.25 * PrimCell.a(2);
    for (double x = -1.0; x <= 1.0; x += 0.0000500113412)
    {
      // for (double x=-0.003; x<=0.003; x+=0.0000011329343481381) {
      P.R[0]    = x * (PrimCell.a(0) + 0.914 * PrimCell.a(1) + 0.781413 * PrimCell.a(2));
      double r  = std::sqrt(dot(P.R[0], P.R[0]));
      double rN = std::sqrt(dot(P.R[0] - N, P.R[0] - N));
      OrbitalSet->evaluate(P, 0, phi, grad, lapl);
      // OrbitalSet->evaluate(P, 0, phi);
      fprintf(fout, "%1.12e ", r * x / std::abs(x));
      for (int j = 0; j < numOrbs; j++)
      {
        double gmag = std::sqrt(dot(grad[j], grad[j]));
        fprintf(fout, "%16.12e ",
                /*phi[j]*phi[j]**/ (-5.0 / r - 0.5 * lapl[j] / phi[j]));
        // double E = -5.0/r -0.5*lapl[j]/phi[j];
        fprintf(fout, "%16.12e ", phi[j]);
        fprintf(fout, "%16.12e ", gmag);
      }
      fprintf(fout, "\n");
    }
    fclose(fout);
  }
#endif
#endif
  //if (sourceName.size() && (ParticleSets.find(sourceName) == ParticleSets.end()))
  //{
  //  app_log() << "  EinsplineSetBuilder creates a ParticleSet " << sourceName << std::endl;
  //  ParticleSet* ions=new ParticleSet;
  //  ions->Lattice=TargetPtcl.Lattice;
  //  ESHDFIonsParser ap(*ions,H5FileID,myComm);
  //  ap.put(XMLRoot);
  //  ap.expand(TileMatrix);
  //  ions->setName(sourceName);
  //  ParticleSets[sourceName]=ions;
  //  //overwrite the lattice and assign random
  //  if(TargetPtcl.Lattice.SuperCellEnum)
  //  {
  //    TargetPtcl.Lattice=ions->Lattice;
  //    makeUniformRandom(TargetPtcl.R);
  //    TargetPtcl.R.setUnit(PosUnit::LatticeUnit);
  //    TargetPtcl.convert2Cart(TargetPtcl.R);
  //    TargetPtcl.createSK();
  //  }
  //}
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
