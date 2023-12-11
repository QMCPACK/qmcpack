//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2023 QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Paul R. C. Kent, kentpr@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


/** @file BsplineReader.cpp
 *
 * Implement super function
 */
#include "EinsplineSetBuilder.h"
#include "BsplineReader.h"
#include "OhmmsData/AttributeSet.h"
#include "Message/CommOperators.h"

#include <array>
#include <filesystem>

namespace qmcplusplus
{
BsplineReader::BsplineReader(EinsplineSetBuilder* e)
    : mybuilder(e), checkNorm(true), saveSplineCoefs(false), rotate(true)
{
  myComm = mybuilder->getCommunicator();
}

BsplineReader::~BsplineReader() = default;

inline std::string make_bandinfo_filename(const std::string& root,
                                          int spin,
                                          int twist,
                                          const Tensor<int, 3>& tilematrix,
                                          int gid)
{
  std::ostringstream oo;
  oo << root << ".tile_" << tilematrix(0, 0) << tilematrix(0, 1) << tilematrix(0, 2) << tilematrix(1, 0)
     << tilematrix(1, 1) << tilematrix(1, 2) << tilematrix(2, 0) << tilematrix(2, 1) << tilematrix(2, 2) << ".spin_"
     << spin << ".tw_" << twist;
  if (gid >= 0)
    oo << ".g" << gid;
  return oo.str();
}


inline std::string make_bandgroup_name(const std::string& root,
                                       int spin,
                                       int twist,
                                       const Tensor<int, 3>& tilematrix,
                                       int first,
                                       int last)
{
  std::ostringstream oo;
  oo << root << ".tile_" << tilematrix(0, 0) << tilematrix(0, 1) << tilematrix(0, 2) << tilematrix(1, 0)
     << tilematrix(1, 1) << tilematrix(1, 2) << tilematrix(2, 0) << tilematrix(2, 1) << tilematrix(2, 2) << ".spin_"
     << spin << ".tw_" << twist << ".l" << first << "u" << last;
  return oo.str();
}

void BsplineReader::setCommon(xmlNodePtr cur)
{
  // check orbital normalization by default
  std::string checkOrbNorm("yes");
  std::string saveCoefs("no");
  OhmmsAttributeSet a;
  a.add(checkOrbNorm, "check_orb_norm");
  a.add(saveCoefs, "save_coefs");
  a.put(cur);

  // allow user to turn off norm check with a warning
  if (checkOrbNorm == "no")
  {
    app_log() << "WARNING: disable orbital normalization check!" << std::endl;
    checkNorm = false;
  }
  saveSplineCoefs = saveCoefs == "yes";
}

std::unique_ptr<SPOSet> BsplineReader::create_spline_set(int spin, xmlNodePtr cur)
{
  int ns(0);
  std::string spo_object_name;
  OhmmsAttributeSet a;
  a.add(ns, "size");
  a.add(spo_object_name, "name");
  a.add(spo_object_name, "id");
  a.put(cur);

  if (ns == 0)
    APP_ABORT_TRACE(__FILE__, __LINE__, "parameter/@size missing");

  if (spo2band.empty())
    spo2band.resize(mybuilder->states.size());

  std::vector<BandInfo>& fullband = (*(mybuilder->FullBands[spin]));

  if (spo2band[spin].empty())
  {
    spo2band[spin].reserve(fullband.size());
    if (!mybuilder->states[spin])
      mybuilder->states[spin] = std::make_unique<SPOSetInfo>();
    mybuilder->clear_states(spin);
    initialize_spo2band(spin, fullband, *mybuilder->states[spin], spo2band[spin]);
  }

  BandInfoGroup vals;
  vals.TwistIndex = fullband[0].TwistIndex;
  vals.GroupID    = 0;
  vals.myName = make_bandgroup_name(mybuilder->getName(), spin, mybuilder->twist_num_, mybuilder->TileMatrix, 0, ns);
  vals.selectBands(fullband, 0, ns, false);

  return create_spline_set(spo_object_name, spin, vals);
}

std::unique_ptr<SPOSet> BsplineReader::create_spline_set(int spin, xmlNodePtr cur, SPOSetInputInfo& input_info)
{
  std::string spo_object_name;
  OhmmsAttributeSet a;
  a.add(spo_object_name, "name");
  a.add(spo_object_name, "id");
  a.put(cur);

  if (spo2band.empty())
    spo2band.resize(mybuilder->states.size());

  std::vector<BandInfo>& fullband = (*(mybuilder->FullBands[spin]));

  if (spo2band[spin].empty())
  {
    spo2band[spin].reserve(fullband.size());
    if (!mybuilder->states[spin])
      mybuilder->states[spin] = std::make_unique<SPOSetInfo>();
    mybuilder->clear_states(spin);
    initialize_spo2band(spin, fullband, *mybuilder->states[spin], spo2band[spin]);
  }

  BandInfoGroup vals;
  vals.TwistIndex = fullband[0].TwistIndex;
  vals.GroupID    = 0;
  vals.myName     = make_bandgroup_name(mybuilder->getName(), spin, mybuilder->twist_num_, mybuilder->TileMatrix,
                                        input_info.min_index(), input_info.max_index());
  vals.selectBands(fullband, spo2band[spin][input_info.min_index()], input_info.max_index() - input_info.min_index(),
                   false);

  return create_spline_set(spo_object_name, spin, vals);
}

/** build index tables to map a state to band with k-point folidng
   * @param bigspace full BandInfo constructed by EinsplineSetBuilder
   * @param sposet SPOSetInfo owned by someone, most likely EinsplinseSetBuilder
   * @param spo2band spo2band[i] is the index in bigspace
   *
   * At gamma or arbitrary kpoints with complex wavefunctions, spo2band[i]==i
   */
void BsplineReader::initialize_spo2band(int spin,
                                        const std::vector<BandInfo>& bigspace,
                                        SPOSetInfo& sposet,
                                        std::vector<int>& spo2band)
{
  spo2band.reserve(bigspace.size());
  int ns = 0;
  for (int i = 0; i < bigspace.size(); ++i)
  {
    spo2band.push_back(i);
    SPOInfo a(ns, bigspace[i].Energy);
    sposet.add(a);
    ns++;
    if (bigspace[i].MakeTwoCopies)
    {
      spo2band.push_back(i);
      SPOInfo b(ns, bigspace[i].Energy);
      sposet.add(b);
      ns++;
    }
  }

  //write to a file
  const Communicate* comm = myComm;
  if (comm->rank())
    return;

  std::filesystem::path aname = make_bandinfo_filename(mybuilder->getName(), spin, mybuilder->twist_num_,
                                                       mybuilder->TileMatrix, comm->getGroupID());
  aname += ".bandinfo.dat";

  std::ofstream o(aname.c_str());
  std::array<char, 1024> s;
  ns            = 0;
  using PosType = QMCTraits::PosType;
  o << "#  Band    State   TwistIndex BandIndex Energy      Kx      Ky      Kz      K1      K2      K3    KmK "
    << std::endl;
  for (int i = 0; i < bigspace.size(); ++i)
  {
    int ti     = bigspace[i].TwistIndex;
    int bi     = bigspace[i].BandIndex;
    double e   = bigspace[i].Energy;
    int nd     = (bigspace[i].MakeTwoCopies) ? 2 : 1;
    PosType k  = mybuilder->PrimCell.k_cart(mybuilder->primcell_kpoints[ti]);
    int s_size = std::snprintf(s.data(), s.size(), "%8d %8d %8d %8d %12.6f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %6d\n",
                               i, ns, ti, bi, e, k[0], k[1], k[2], mybuilder->primcell_kpoints[ti][0],
                               mybuilder->primcell_kpoints[ti][1], mybuilder->primcell_kpoints[ti][2], nd);
    if (s_size < 0)
      throw std::runtime_error("Error generating bandinfo");
    o << s.data();
    ns += nd;
  }
}
} // namespace qmcplusplus
