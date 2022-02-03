//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2017 Jeongnim Kim and QMCPACK developers.
//
// File developed by:  Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#undef NDEBUG

#include "catch.hpp"

#include "Configuration.h"

// Avoid the need to link with other libraries just to get APP_ABORT
#undef APP_ABORT
#define APP_ABORT(x)             \
  {                              \
    std::cout << x << std::endl; \
    throw;                       \
  }

#include "OhmmsData/Libxml2Doc.h"
#include "Utilities/RandomGenerator.h"
#include "ProjectData.h"

#include "hdf/hdf_multi.h"
#include "hdf/hdf_archive.h"

#include <stdio.h>
#include <string>
#include <vector>
#include <complex>

#include "mpi3/communicator.hpp"
#include "mpi3/shared_communicator.hpp"
//#include "mpi3/environment.hpp"

//#include "AFQMC/Walkers WalkerSetFactory.hpp"
#include "AFQMC/Walkers/WalkerSet.hpp"
#include "AFQMC/Walkers/WalkerIO.hpp"

using std::complex;
using std::cout;
using std::endl;
using std::string;

namespace qmcplusplus
{
using namespace afqmc;

void myREQUIRE(const double& a, const double& b) { REQUIRE(a == Approx(b)); }

void myREQUIRE(const std::complex<double>& a, const double& b) { REQUIRE(a.real() == Approx(b)); }

void myREQUIRE(const std::complex<double>& a, const std::complex<double>& b)
{
  REQUIRE(a.real() == Approx(b.real()));
  REQUIRE(a.imag() == Approx(b.imag()));
}

template<class M1, class M2>
void check(M1&& A, M2& B)
{
  using element1 = typename std::decay<M1>::type::element;
  using element2 = typename std::decay<M2>::type::element;
  REQUIRE(A.size(0) == B.size(0));
  REQUIRE(A.size(1) == B.size(1));
  for (int i = 0; i < A.size(0); i++)
    for (int j = 0; j < A.size(1); j++)
      myREQUIRE(element1(A[i][j]), element2(B[i][j]));
}

using namespace afqmc;
using communicator = boost::mpi3::communicator;

void test_basic_walker_features(bool serial, std::string wtype)
{
  auto world = boost::mpi3::environment::get_world_instance();
  auto node  = world.split_shared(world.rank());

#if defined(ENABLE_CUDA) || defined(ENABLE_HIP)
  arch::INIT(node);
#endif

  using Type = std::complex<double>;

  //assert(world.size()%2 == 0);

  int NMO = 8, NAEA = 2, NAEB = 2, nwalkers = 10;
  if (wtype == "noncollinear")
  {
    NAEA = 4;
    NAEB = 0;
  }

  //auto node = world.split_shared();

  GlobalTaskGroup gTG(world);
  TaskGroup_ TG(gTG, std::string("TaskGroup"), 1, serial ? 1 : gTG.getTotalCores());
  AFQMCInfo info;
  info.NMO  = NMO;
  info.NAEA = NAEA;
  info.NAEB = NAEB;
  info.name = "walker";
  int M((wtype == "noncollinear") ? 2 * NMO : NMO);
  boost::multi::array<Type, 2> initA({M, NAEA});
  boost::multi::array<Type, 2> initB({M, NAEB});
  for (int i = 0; i < NAEA; i++)
    initA[i][i] = Type(0.22);
  for (int i = 0; i < NAEB; i++)
    initB[i][i] = Type(0.22);
  RandomGenerator rng;

  std::string xml_block;
  xml_block = "<WalkerSet name=\"wset0\">  \
  <parameter name=\"min_weight\">0.05</parameter>  \
  <parameter name=\"max_weight\">4</parameter>  \
  <parameter name=\"walker_type\">" +
      wtype + "</parameter>  \
  <parameter name=\"load_balance\">async</parameter>  \
  <parameter name=\"pop_control\">pair</parameter>  \
</WalkerSet> \
";
  Libxml2Document doc;
  bool okay = doc.parseFromString(xml_block.c_str());
  REQUIRE(okay);

  WalkerSet wset(TG, doc.getRoot(), info, &rng);
  wset.resize(nwalkers, initA, initB);

  REQUIRE(wset.size() == nwalkers);
  int cnt           = 0;
  double tot_weight = 0.0;
  for (WalkerSet::iterator it = wset.begin(); it != wset.end(); ++it)
  {
    auto sm = it->SlaterMatrix(Alpha);
    REQUIRE(*it->SlaterMatrix(Alpha) == initA);
    *it->weight()  = cnt * 1.0 + 0.5;
    *it->overlap() = cnt * 1.0 + 0.5;
    *it->E1()      = cnt * 1.0 + 0.5;
    *it->EXX()     = cnt * 1.0 + 0.5;
    *it->EJ()      = cnt * 1.0 + 0.5;
    tot_weight += cnt * 1.0 + 0.5;
    cnt++;
  }
  REQUIRE(cnt == nwalkers);
  cnt = 0;
  for (WalkerSet::iterator it = wset.begin(); it != wset.end(); ++it)
  {
    Type d_(cnt * 1.0 + 0.5);
    REQUIRE(*it->weight() == d_);
    REQUIRE(*it->overlap() == cnt * 1.0 + 0.5);
    REQUIRE(*it->E1() == cnt * 1.0 + 0.5);
    REQUIRE(*it->EXX() == cnt * 1.0 + 0.5);
    REQUIRE(*it->EJ() == cnt * 1.0 + 0.5);
    cnt++;
  }

  wset.reserve(20);
  REQUIRE(wset.capacity() == 20);
  cnt = 0;
  for (WalkerSet::iterator it = wset.begin(); it != wset.end(); ++it)
  {
    REQUIRE(*it->weight() == cnt * 1.0 + 0.5);
    REQUIRE(*it->overlap() == cnt * 1.0 + 0.5);
    REQUIRE(*it->E1() == cnt * 1.0 + 0.5);
    REQUIRE(*it->EXX() == cnt * 1.0 + 0.5);
    REQUIRE(*it->EJ() == cnt * 1.0 + 0.5);
    cnt++;
  }
  for (int i = 0; i < wset.size(); i++)
  {
    REQUIRE(*wset[i].weight() == i * 1.0 + 0.5);
    REQUIRE(*wset[i].overlap() == i * 1.0 + 0.5);
    REQUIRE(*wset[i].E1() == i * 1.0 + 0.5);
    REQUIRE(*wset[i].EXX() == i * 1.0 + 0.5);
    REQUIRE(*wset[i].EJ() == i * 1.0 + 0.5);
  }
  for (int i = 0; i < wset.size(); i++)
  {
    auto w = wset[i];
    REQUIRE(*w.weight() == i * 1.0 + 0.5);
    REQUIRE(*w.overlap() == i * 1.0 + 0.5);
    REQUIRE(*w.E1() == i * 1.0 + 0.5);
    REQUIRE(*w.EXX() == i * 1.0 + 0.5);
    REQUIRE(*w.EJ() == i * 1.0 + 0.5);
  }
  REQUIRE(wset.get_TG_target_population() == nwalkers);
  REQUIRE(wset.get_global_target_population() == nwalkers * TG.getNumberOfTGs());
  REQUIRE(wset.GlobalPopulation() == nwalkers * TG.getNumberOfTGs());
  REQUIRE(wset.GlobalPopulation() == wset.get_global_target_population());
  REQUIRE(wset.NumBackProp() == 0);
  REQUIRE(wset.GlobalWeight() == tot_weight * TG.getNumberOfTGs());

  wset.scaleWeight(2.0);
  tot_weight *= 2.0;
  REQUIRE(wset.GlobalWeight() == tot_weight * TG.getNumberOfTGs());

  std::vector<ComplexType> Wdata;
  wset.popControl(Wdata);
  REQUIRE(wset.GlobalWeight() == Approx(static_cast<RealType>(wset.get_global_target_population())));
  REQUIRE(wset.get_TG_target_population() == nwalkers);
  REQUIRE(wset.get_global_target_population() == nwalkers * TG.getNumberOfTGs());
  REQUIRE(wset.GlobalPopulation() == nwalkers * TG.getNumberOfTGs());
  REQUIRE(wset.GlobalPopulation() == wset.get_global_target_population());
  REQUIRE(wset.GlobalWeight() == Approx(static_cast<RealType>(wset.get_global_target_population())));
  double nx = (wset.getWalkerType() == NONCOLLINEAR ? 1.0 : 2.0);
  for (int i = 0; i < wset.size(); i++)
  {
    auto w = wset[i];
    myREQUIRE(std::exp(nx * wset.getLogOverlapFactor()) * ComplexType(*w.overlap()), ComplexType(*w.E1()));
    REQUIRE(*w.EXX() == *w.E1());
    REQUIRE(*w.EJ() == *w.E1());
  }

  wset.clean();
  REQUIRE(wset.size() == 0);
  REQUIRE(wset.capacity() == 0);
}

void test_hyperslab()
{
  auto world = boost::mpi3::environment::get_world_instance();
  auto node  = world.split_shared(world.rank());

#if defined(ENABLE_CUDA) || defined(ENABLE_HIP)
  arch::INIT(node);
#endif

  using Type   = std::complex<double>;
  using Matrix = boost::multi::array<Type, 2>;

  int rank = world.rank();

  int nwalk         = 9;
  int nprop         = 7;
  int nprop_to_safe = 3;
  Matrix Data({nwalk, nprop});

  for (int i = 0; i < nwalk; i++)
    for (int j = 0; j < nprop; j++)
      Data[i][j] = i * 10 + rank * 100 + j;

  int nwtot = (world += nwalk);

  hdf_archive dump(world, true);
  if (!dump.create("dummy_walkers.h5", H5F_ACC_EXCL))
  {
    app_error() << " Error opening restart file. \n";
    APP_ABORT("");
  }
  dump.push("WalkerSet");

  hyperslab_proxy<Matrix, 2> hslab(Data, std::array<size_t, 2>{static_cast<size_t>(nwtot), static_cast<size_t>(nprop)},
                                   std::array<size_t, 2>{static_cast<size_t>(nwalk), static_cast<size_t>(nprop)},
                                   std::array<size_t, 2>{static_cast<size_t>(rank * nwalk), 0});
  dump.write(hslab, "Walkers");
  dump.close();
  world.barrier();

  {
    hdf_archive read(world, false);
    if (!read.open("dummy_walkers.h5", H5F_ACC_RDONLY))
    {
      app_error() << " Error opening restart file. \n";
      APP_ABORT("");
    }
    read.push("WalkerSet");

    Matrix DataIn({nwalk, nprop});

    hyperslab_proxy<Matrix, 2> hslab(DataIn,
                                     std::array<size_t, 2>{static_cast<size_t>(nwtot), static_cast<size_t>(nprop)},
                                     std::array<size_t, 2>{static_cast<size_t>(nwalk), static_cast<size_t>(nprop)},
                                     std::array<size_t, 2>{static_cast<size_t>(rank * nwalk), 0});
    read.read(hslab, "Walkers");
    read.close();

    for (int i = 0; i < nwalk; i++)
      for (int j = 0; j < nprop; j++)
      {
        REQUIRE(real(DataIn[i][j]) == i * 10 + rank * 100 + j);
        REQUIRE(imag(DataIn[i][j]) == 0);
      }
  }
  world.barrier();
  if (world.root())
    remove("dummy_walkers.h5");
}

void test_double_hyperslab()
{
  auto world = boost::mpi3::environment::get_world_instance();

  using Type   = std::complex<double>;
  using Matrix = boost::multi::array<Type, 2>;

  int rank = world.rank();

  int nwalk         = 9;
  int nprop         = 3;
  int nprop_to_safe = 3;
  Matrix Data({nwalk, nprop});

  for (int i = 0; i < nwalk; i++)
    for (int j = 0; j < nprop; j++)
      Data[i][j] = i * 10 + rank * 100 + j;

  int nwtot = (world += nwalk);

  hdf_archive dump(world, true);
  if (!dump.create("dummy_walkers.h5", H5F_ACC_EXCL))
  {
    app_error() << " Error opening restart file. \n";
    APP_ABORT("");
  }
  dump.push("WalkerSet");

  //double_hyperslab_proxy<Matrix,2> hslab(Data,
  hyperslab_proxy<Matrix, 2> hslab(Data,
                                   std::array<size_t, 2>{static_cast<size_t>(nwtot),
                                                         static_cast<size_t>(nprop_to_safe)},
                                   std::array<size_t, 2>{static_cast<size_t>(nwalk),
                                                         static_cast<size_t>(nprop_to_safe)},
                                   std::array<size_t, 2>{static_cast<size_t>(rank * nwalk), 0}); //,

  //                                  std::array<int,2>{nwalk,nprop},
  //                                  std::array<int,2>{nwalk,nprop_to_safe},
  //                                  std::array<int,2>{0,0});
  dump.write(hslab, "Walkers");
  dump.close();
  world.barrier();

  {
    hdf_archive read(world, false);
    if (!read.open("dummy_walkers.h5", H5F_ACC_RDONLY))
    {
      app_error() << " Error opening restart file. \n";
      APP_ABORT("");
    }
    read.push("WalkerSet");

    //Matrix DataIn({nwalk,nprop});
    Matrix DataIn({nwalk, nprop_to_safe});

    //double_hyperslab_proxy<Matrix,2> hslab(DataIn,
    hyperslab_proxy<Matrix, 2> hslab(DataIn,
                                     std::array<size_t, 2>{static_cast<size_t>(nwtot),
                                                           static_cast<size_t>(nprop_to_safe)},
                                     std::array<size_t, 2>{static_cast<size_t>(nwalk),
                                                           static_cast<size_t>(nprop_to_safe)},
                                     std::array<size_t, 2>{static_cast<size_t>(rank * nwalk), 0}); //,
    //                                  std::array<int,2>{nwalk,nprop},
    //                                  std::array<int,2>{nwalk,nprop_to_safe},
    //                                  std::array<int,2>{0,0});
    read.read(hslab, "Walkers");
    read.close();

    for (int i = 0; i < nwalk; i++)
    {
      for (int j = 0; j < nprop_to_safe; j++)
      {
        REQUIRE(real(DataIn[i][j]) == i * 10 + rank * 100 + j);
        REQUIRE(imag(DataIn[i][j]) == 0);
      }
      /*
     for(int j=nprop_to_safe; j<nprop; j++) {
       REQUIRE( real(DataIn[i][j]) == 0);
       REQUIRE( imag(DataIn[i][j]) == 0);
     }
*/
    }
  }
  world.barrier();
  if (world.root())
    remove("dummy_walkers.h5");
}

void test_walker_io(std::string wtype)
{
  auto world = boost::mpi3::environment::get_world_instance();
  auto node  = world.split_shared(world.rank());

  using Type = std::complex<double>;

#if defined(ENABLE_CUDA) || defined(ENABLE_HIP)
  arch::INIT(node);
#endif

  //assert(world.size()%2 == 0);

  int NMO = 8, NAEA = 2, NAEB = 2, nwalkers = 10;
  if (wtype == "noncollinear")
  {
    NAEA = 4;
    NAEB = 0;
  }

  //auto node = world.split_shared();

  GlobalTaskGroup gTG(world);
  TaskGroup_ TG(gTG, std::string("TaskGroup"), 1, 1);
  AFQMCInfo info;
  info.NMO  = NMO;
  info.NAEA = NAEA;
  info.NAEB = NAEB;
  info.name = "walker";
  int M((wtype == "noncollinear") ? 2 * NMO : NMO);
  boost::multi::array<Type, 2> initA({M, NAEA});
  boost::multi::array<Type, 2> initB({M, NAEB});
  for (int i = 0; i < NAEA; i++)
    initA[i][i] = Type(0.22);
  for (int i = 0; i < NAEB; i++)
    initB[i][i] = Type(0.22);
  RandomGenerator rng;

  std::string xml_block;
  xml_block = "<WalkerSet name=\"wset0\">  \
  <parameter name=\"walker_type\">" +
      wtype + "</parameter>  \
</WalkerSet> \
";
  Libxml2Document doc;
  bool okay = doc.parseFromString(xml_block.c_str());
  REQUIRE(okay);

  WalkerSet wset(TG, doc.getRoot(), info, &rng);
  wset.resize(nwalkers, initA, initB);

  REQUIRE(wset.size() == nwalkers);
  int cnt = 0;
  for (WalkerSet::iterator it = wset.begin(); it != wset.end(); ++it)
  {
    auto sm = it->SlaterMatrix(Alpha);
    REQUIRE(*it->SlaterMatrix(Alpha) == initA);
    *it->weight()  = cnt * 1.0 + 0.5;
    *it->overlap() = cnt * 1.0 + 0.5;
    *it->E1()      = cnt * 1.0 + 0.5;
    *it->EXX()     = cnt * 1.0 + 0.5;
    *it->EJ()      = cnt * 1.0 + 0.5;
    cnt++;
  }
  REQUIRE(cnt == nwalkers);

#if defined(ENABLE_PHDF5)
  hdf_archive dump(world, true);
  {
#else
  hdf_archive dump(world, false);
  if (TG.Global().root())
  {
#endif
    if (!dump.create("dummy_walkers.h5", H5F_ACC_EXCL))
    {
      app_error() << " Error opening restart file. \n";
      APP_ABORT("");
    }
  }

  // dump restart file
  dumpToHDF5(wset, dump);
  dump.close();

  {
#if defined(ENABLE_PHDF5)
    hdf_archive read(world, true);
    {
#else
    hdf_archive read(world, false);
    if (TG.Global().root())
    {
#endif
      if (!read.open("dummy_walkers.h5", H5F_ACC_RDONLY))
      {
        app_error() << " Error opening restart file. \n";
        APP_ABORT("");
      }
      else
      {
        read.close();
      }
    }

    WalkerSet wset2(TG, doc.getRoot(), info, &rng);
    restartFromHDF5(wset2, nwalkers, "dummy_walkers.h5", read, true);
    for (int i = 0; i < nwalkers; i++)
    {
      REQUIRE(*wset[i].SlaterMatrix(Alpha) == *wset2[i].SlaterMatrix(Alpha));
      REQUIRE(*wset[i].weight() == *wset2[i].weight());
      REQUIRE(*wset[i].overlap() == *wset2[i].overlap());
      REQUIRE(*wset[i].E1() == *wset2[i].E1());
      REQUIRE(*wset[i].EXX() == *wset2[i].EXX());
      REQUIRE(*wset[i].EJ() == *wset2[i].EJ());
    }
  }
  world.barrier();
  if (world.root())
    remove("dummy_walkers.h5");
}

TEST_CASE("swset_test_serial", "[shared_wset]")
{
  test_basic_walker_features(true, "closed");
  test_basic_walker_features(false, "closed");
  test_basic_walker_features(true, "collinear");
  test_basic_walker_features(false, "collinear");
  test_basic_walker_features(true, "noncollinear");
  test_basic_walker_features(false, "noncollinear");
}
/*
TEST_CASE("hyperslab_tests", "[shared_wset]")
{
 // test_hyperslab();
  test_double_hyperslab();
}
*/
TEST_CASE("walker_io", "[shared_wset]")
{
  test_walker_io("closed");
  test_walker_io("collinear");
  test_walker_io("noncollinear");
}

} // namespace qmcplusplus
