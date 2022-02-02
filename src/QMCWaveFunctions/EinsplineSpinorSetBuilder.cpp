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


#include "EinsplineSpinorSetBuilder.h"
#include "QMCWaveFunctions/SpinorSet.h"
#include "OhmmsData/AttributeSet.h"
#include "Message/CommOperators.h"
#include "Utilities/Timer.h"
#include "QMCWaveFunctions/einspline_helper.hpp"
#include "QMCWaveFunctions/BsplineFactory/BsplineReaderBase.h"
#include "QMCWaveFunctions/BsplineFactory/createBsplineReader.h"

namespace qmcplusplus
{
std::unique_ptr<SPOSet> EinsplineSpinorSetBuilder::createSPOSetFromXML(xmlNodePtr cur)
{
  int numOrbs = 0;
  int sortBands(1);
  int spinSet       = 0;
  int spinSet2      = 1;
  int twist_num_inp = TWISTNUM_NO_INPUT;
  TinyVector<double, OHMMS_DIM> twist_inp(TWIST_NO_INPUT);

  //There have to be two "spin states"...  one for the up channel and one for the down channel.
  // We force this for spinors and manually resize states and FullBands.
  states.clear();
  states.resize(2);

  FullBands.resize(2);

  SPOSet* UpOrbitalSet;
  std::string sourceName;
  std::string spo_prec("double");
  std::string truncate("no");
  std::string hybrid_rep("no");
  ScopedTimer spo_timer_scope(*timer_manager.createTimer("einspline::CreateSpinorSetFromXML", timer_level_medium));

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
    a.add(spo_prec, "precision");
    a.add(truncate, "truncate");
    a.add(myName, "tag");

    a.put(XMLRoot);
    a.add(numOrbs, "size");
    a.add(numOrbs, "norbs");
    a.add(spinSet, "spindataset");
    a.add(spinSet, "group");
    a.put(cur);
  }

  auto pit(ParticleSets.find(sourceName));
  if (pit == ParticleSets.end())
    myComm->barrier_and_abort("Einspline needs the source particleset");
  else
    SourcePtcl = pit->second;

  ///////////////////////////////////////////////
  // Read occupation information from XML file //
  ///////////////////////////////////////////////
  std::vector<int> Occ_Old(0, 0);
  Occ.resize(0, 0);
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
        myComm->barrier_and_abort("EinsplineSetBuilder::createSPOSet Only ground state occupation currently "
                                  "supported in EinsplineSetBuilder.");
    }
    cur = cur->next;
  }
  if (Occ != Occ_Old)
  {
    NewOcc  = true;
    Occ_Old = Occ;
  }
  else
    NewOcc = false;

  H5OrbSet aset(H5FileName, spinSet, numOrbs);
  std::map<H5OrbSet, SPOSet*, H5OrbSet>::iterator iter;
  iter = SPOSetMap.find(aset);
  if ((iter != SPOSetMap.end()) && (!NewOcc))
  {
    app_log() << "SPOSet parameters match in EinsplineSetBuilder:  "
              << "cloning EinsplineSet object.\n";
    return std::unique_ptr<SPOSet>{iter->second->makeClone()};
  }

  if (FullBands[spinSet] == 0)
    FullBands[spinSet] = std::make_unique<std::vector<BandInfo>>();

  if (FullBands[spinSet2] == 0)
    FullBands[spinSet2] = std::make_unique<std::vector<BandInfo>>();

  //This is to skip checks on ion-ID's, spin types, etc.  If we've made it here, we assume we know better
  //than Einspline on what the data means...
  bool skipChecks = true;

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

  std::string useGPU("no");
#if !defined(QMC_COMPLEX)
  if (use_real_splines_)
  {
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
  //Norm for spinor wavefunctions is different from SPO's by a factor of sqrt(2).  Disable the unit norm check.
  MixedSplineReader->setCheckNorm(false);
  //Set no rotation to the orbitals
  MixedSplineReader->setRotate(false);

  //Make the up spin set.
  HasCoreOrbs       = bcastSortBands(spinSet, NumDistinctOrbitals, myComm->rank() == 0);
  auto bspline_zd_u = MixedSplineReader->create_spline_set(spinSet, spo_cur);

  //Make the down spin set.
  OccupyBands(spinSet2, sortBands, numOrbs, skipChecks);
  HasCoreOrbs       = bcastSortBands(spinSet2, NumDistinctOrbitals, myComm->rank() == 0);
  auto bspline_zd_d = MixedSplineReader->create_spline_set(spinSet2, spo_cur);

  //register with spin set and we're off to the races.
  auto spinor_set = std::make_unique<SpinorSet>();
  spinor_set->set_spos(std::move(bspline_zd_u), std::move(bspline_zd_d));
  return spinor_set;
};
} // namespace qmcplusplus
