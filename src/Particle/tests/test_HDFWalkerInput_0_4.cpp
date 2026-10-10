//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2026 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#include <catch2/catch_test_macros.hpp>
#include "Utilities/for_testing/Catch2Approx.h"

#include <filesystem>
#include <string>
#include <vector>

#include "Message/Communicate.h"
#include "Message/CommOperators.h"
#include "Particle/WalkerConfigurations.h"
#include "Particle/HDFWalkerOutput.h"
#include "Particle/HDFWalkerInput_0_4.h"
#include "hdf/hdf_archive.h"
#include "hdf/HDFVersion.h"
#include "Utilities/FairDivide.h"
#include "OhmmsData/Libxml2Doc.h"

namespace qmcplusplus
{

namespace
{
constexpr size_t total_walkers = 16;
constexpr size_t num_ptcls     = 2;

inline QMCTraits::RealType expected_R(size_t walker_id, size_t ptcl_id, size_t dim)
{
  return static_cast<QMCTraits::RealType>(walker_id * 100.0 + ptcl_id * 10.0 + dim + 1.0);
}

inline QMCTraits::FullPrecRealType expected_weight(size_t walker_id)
{
  return static_cast<QMCTraits::FullPrecRealType>(1.0 + 0.125 * walker_id);
}

void verify_walkers(const WalkerConfigurations& wc_list,
                    size_t expected_count,
                    size_t global_offset,
                    size_t nptcls)
{
  REQUIRE(wc_list.getActiveWalkers() == expected_count);
  for (size_t iw = 0; iw < expected_count; ++iw)
  {
    const size_t g = global_offset + iw;
    for (size_t ip = 0; ip < nptcls; ++ip)
      for (size_t id = 0; id < OHMMS_DIM; ++id)
        CHECK(wc_list[iw]->R[ip][id] == Approx(expected_R(g, ip, id)));
    CHECK(wc_list[iw]->Weight == Approx(expected_weight(g)));
  }
}
} // namespace

TEST_CASE("HDFWalkerInput_0_4 16 walkers read_hdf5 and put", "[particle]")
{
  Communicate& c(*OHMMS::Controller);

  std::vector<int> walker_offsets(c.size() + 1);
  FairDivideLow(total_walkers, c.size(), walker_offsets);
  const int nw_local = walker_offsets[c.rank() + 1] - walker_offsets[c.rank()];

  // Create local walkers
  WalkerConfigurations wc_list_out;
  wc_list_out.createWalkers(nw_local, num_ptcls);
  for (int iw = 0; iw < nw_local; ++iw)
  {
    const size_t g = walker_offsets[c.rank()] + iw;
    for (size_t ip = 0; ip < num_ptcls; ++ip)
      for (size_t id = 0; id < OHMMS_DIM; ++id)
        wc_list_out[iw]->R[ip][id] = expected_R(g, ip, id);
    wc_list_out[iw]->Weight = expected_weight(g);
  }
  wc_list_out.setWalkerOffsets(walker_offsets);

  const std::string file_root = "test_hdfwalkerinput_16";
  c.setName(file_root);
  HDFWalkerOutput hout(num_ptcls, file_root, c);
  hout.dump(wc_list_out, 0);
  c.barrier();

  const std::filesystem::path h5file = file_root + ".config.h5";
  REQUIRE(std::filesystem::exists(h5file));

  const HDFVersion version(0, 4);

  // Test 1: read_hdf5
  {
    WalkerConfigurations wc_list_in;
    HDFWalkerInput_0_4 hinp(wc_list_in, num_ptcls, c, version);
    const bool success = hinp.read_hdf5(h5file);
    REQUIRE(success);

    int total_read = wc_list_in.getActiveWalkers();
    c.allreduce(total_read);
    CHECK(total_read == total_walkers);

    verify_walkers(wc_list_in, nw_local, walker_offsets[c.rank()], num_ptcls);
  }

  // Test 2: put via XML
  {
    const std::string xml_str = "<mcwalkerset fileroot=\"" + file_root + "\" version=\"0 4\" collected=\"yes\"/>";
    Libxml2Document doc;
    const bool parsed = doc.parseFromString(xml_str);
    REQUIRE(parsed);

    WalkerConfigurations wc_list_xml;
    HDFWalkerInput_0_4 hinp_xml(wc_list_xml, num_ptcls, c, version);
    const bool success = hinp_xml.put(doc.getRoot());
    REQUIRE(success);

    int total_read = wc_list_xml.getActiveWalkers();
    c.allreduce(total_read);
    CHECK(total_read == total_walkers);

    verify_walkers(wc_list_xml, nw_local, walker_offsets[c.rank()], num_ptcls);
  }

  // Test 3: read_hdf5_scatter
  {
    WalkerConfigurations wc_list_scatter;
    HDFWalkerInput_0_4 hinp_scatter(wc_list_scatter, num_ptcls, c, version);
    const bool success = hinp_scatter.read_hdf5_scatter(h5file);
    REQUIRE(success);

    int total_read = wc_list_scatter.getActiveWalkers();
    c.allreduce(total_read);
    CHECK(total_read == total_walkers);

    verify_walkers(wc_list_scatter, nw_local, walker_offsets[c.rank()], num_ptcls);
  }

#if defined(ENABLE_PHDF5)
  // Test 4: read_phdf5
  {
    WalkerConfigurations wc_list_phdf5;
    HDFWalkerInput_0_4 hinp_phdf5(wc_list_phdf5, num_ptcls, c, version);
    const bool success = hinp_phdf5.read_phdf5(h5file);
    REQUIRE(success);

    int total_read = wc_list_phdf5.getActiveWalkers();
    c.allreduce(total_read);
    CHECK(total_read == total_walkers);

    verify_walkers(wc_list_phdf5, nw_local, walker_offsets[c.rank()], num_ptcls);
  }
#endif

  c.barrier();
}

TEST_CASE("HDFWalkerInput_0_4 16 walkers repartition", "[particle]")
{
  Communicate& c(*OHMMS::Controller);

  const std::string single_root = "test_hdfwalkerinput_16_single";
  const std::filesystem::path h5file = single_root + ".config.h5";

  // Master rank writes 16 walkers with a single partition [0, 16]
  if (c.rank() == 0)
  {
    hdf_archive hout(c, false);
    hout.create(h5file);
    HDFVersion cur_version(0, 4);
    hout.write(cur_version.version, hdf::version);
    hout.push(hdf::main_state);
    hout.write(0, "block");
    const size_t nw = total_walkers;
    hout.write(nw, hdf::num_walkers);

    std::vector<int> single_partition{0, static_cast<int>(total_walkers)};
    hout.write(single_partition, "walker_partition");

    std::vector<QMCTraits::RealType> pos_data(total_walkers * num_ptcls * OHMMS_DIM);
    std::vector<QMCTraits::FullPrecRealType> wt_data(total_walkers);
    for (size_t g = 0; g < total_walkers; ++g)
    {
      wt_data[g] = expected_weight(g);
      for (size_t ip = 0; ip < num_ptcls; ++ip)
        for (size_t id = 0; id < OHMMS_DIM; ++id)
          pos_data[g * num_ptcls * OHMMS_DIM + ip * OHMMS_DIM + id] = expected_R(g, ip, id);
    }
    std::array<size_t, 3> gcounts{total_walkers, num_ptcls, OHMMS_DIM};
    hout.writeSlabReshaped(pos_data, gcounts, hdf::walkers);
    std::array<size_t, 1> wcounts{total_walkers};
    hout.writeSlabReshaped(wt_data, wcounts, hdf::walker_weights);
    hout.close();
  }
  c.barrier();

  std::vector<int> expected_offsets(c.size() + 1);
  FairDivideLow(total_walkers, c.size(), expected_offsets);
  const int nw_local = expected_offsets[c.rank() + 1] - expected_offsets[c.rank()];

  const HDFVersion version(0, 4);

  // read_hdf5 repartitioning
  {
    WalkerConfigurations wc_list_repart;
    HDFWalkerInput_0_4 hinp(wc_list_repart, num_ptcls, c, version);
    const bool success = hinp.read_hdf5(h5file);
    REQUIRE(success);

    int total_read = wc_list_repart.getActiveWalkers();
    c.allreduce(total_read);
    CHECK(total_read == total_walkers);

    verify_walkers(wc_list_repart, nw_local, expected_offsets[c.rank()], num_ptcls);
  }

  // read_hdf5_scatter repartitioning
  {
    WalkerConfigurations wc_list_scatter;
    HDFWalkerInput_0_4 hinp_scatter(wc_list_scatter, num_ptcls, c, version);
    const bool success = hinp_scatter.read_hdf5_scatter(h5file);
    REQUIRE(success);

    int total_read = wc_list_scatter.getActiveWalkers();
    c.allreduce(total_read);
    CHECK(total_read == total_walkers);

    verify_walkers(wc_list_scatter, nw_local, expected_offsets[c.rank()], num_ptcls);
  }

#if defined(ENABLE_PHDF5)
  // read_phdf5 repartitioning
  {
    WalkerConfigurations wc_list_phdf5;
    HDFWalkerInput_0_4 hinp_phdf5(wc_list_phdf5, num_ptcls, c, version);
    const bool success = hinp_phdf5.read_phdf5(h5file);
    REQUIRE(success);

    int total_read = wc_list_phdf5.getActiveWalkers();
    c.allreduce(total_read);
    CHECK(total_read == total_walkers);

    verify_walkers(wc_list_phdf5, nw_local, expected_offsets[c.rank()], num_ptcls);
  }
#endif

  c.barrier();
}

TEST_CASE("HDFWalkerInput_0_4 error handling", "[particle]")
{
  Communicate& c(*OHMMS::Controller);

  const HDFVersion version(0, 4);
  WalkerConfigurations wc_list;
  HDFWalkerInput_0_4 hinp(wc_list, num_ptcls, c, version);

  // Reading non-existent file returns false
  const bool success_nonexistent = hinp.read_hdf5("non_existent_file.config.h5");
  CHECK(!success_nonexistent);

  const bool success_scatter_nonexistent = hinp.read_hdf5_scatter("non_existent_file.config.h5");
  CHECK(!success_scatter_nonexistent);

  // Put with empty or invalid node returns false
  const std::string xml_str = "<mcwalkerset/>";
  Libxml2Document doc;
  const bool parsed = doc.parseFromString(xml_str);
  REQUIRE(parsed);
  const bool success_empty_xml = hinp.put(doc.getRoot());
  CHECK(!success_empty_xml);

  c.barrier();
}

} // namespace qmcplusplus
