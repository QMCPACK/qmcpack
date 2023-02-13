//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//
// File created by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#undef NDEBUG

#include "catch.hpp"

//#include "catch.hpp"
#include "Configuration.h"

#include "ProjectData.h"

// Avoid the need to link with other libraries just to get APP_ABORT
#undef APP_ABORT
#define APP_ABORT(x)             \
  {                              \
    std::cout << x << std::endl; \
    throw;                       \
  }

#include <stdio.h>
#include <string>
#include <vector>
#include <complex>

#include "mpi3/shared_communicator.hpp"
#include "mpi3/environment.hpp"

#include "multi/array.hpp"
#include "multi/array_ref.hpp"

#include "AFQMC/Memory/buffer_managers.h"
#include "AFQMC/Matrix/csr_matrix_construct.hpp"
#include "AFQMC/SlaterDeterminantOperations/SlaterDetOperations.hpp"
#include "AFQMC/SlaterDeterminantOperations/mixed_density_matrix.hpp"

using std::complex;
using std::cout;
using std::endl;
using std::string;

namespace qmcplusplus
{
namespace afqmc
{
void myCHECK(const double& a, const double& b) { CHECK(a == Approx(b)); }

void myCHECK(const std::complex<double>& a, const std::complex<double>& b)
{
  CHECK(a.real() == Approx(b.real()));
  CHECK(a.imag() == Approx(b.imag()));
}

template<class M1, class M2>
void check(M1&& A, M2& B)
{
  REQUIRE(std::get<0>(A.sizes()) == std::get<0>(B.sizes()));
  REQUIRE(std::get<1>(A.sizes()) == std::get<1>(B.sizes()));
  for (int i = 0; i < std::get<0>(A.sizes()); i++)
    for (int j = 0; j < std::get<1>(A.sizes()); j++)
      myCHECK(A[i][j], B[i][j]);
}

using namespace afqmc;
/*
TEST_CASE("SDetOps_double_serial", "[sdet_ops]")
{
  Communicate *c;
  //c = OHMMS::Controller;

  const int NMO = 4;
  const int NEL = 3;

  using Type = RealType;
  using vector = std::vector<Type>;
  using array = boost::multi::array<Type,2>;
  using array_ref = boost::multi::array_ref<Type,2>;

  const Type ov = 5.10443199999999;
  const Type ov2 = -11.0204000000000;

  // some arbitrary matrices
  // actually need transpose
  vector m_a = {
//   0.90000,   2.40000,   3.00000,
//   0.40000,   1.00000,   1.20000,
//   1.40000,   1.60000,   3.60000,
//   0.40000,   0.20000,   0.10000
    0.90000, 0.40000, 1.40000, 0.40000,
    2.40000, 1.00000, 1.60000, 0.20000,
    3.00000, 1.20000, 3.60000, 0.10000
  };
  vector m_b = {
   1.90000,   1.40000,   0.40000,
   1.40000,   0.20000,   2.20000,
   0.40000,   2.60000,   0.60000,
   1.10000,   0.30000,   0.90000
  };

  array A({NEL,NMO});
  array B({NMO,NEL});

  for(int i=0, k=0; i<A.size(0); i++)
    for(int j=0; j<A.size(1); j++,k++)
       A[i][j] = m_a[k];

  for(int i=0, k=0; i<B.size(0); i++)
    for(int j=0; j<B.size(1); j++,k++)
       B[i][j] = m_b[k];

  array_ref Aref(m_a.data(),{NEL,NMO});
  array_ref Bref(m_b.data(),{NMO,NEL});

  SlaterDetOperations SDet( SlaterDetOperations_shared<Type>(NMO,NEL) );

  // Overlaps
  CHECK(SDet.Overlap(A,B) == Approx(ov));
  CHECK(SDet.Overlap(Aref,B) == Approx(ov));
  CHECK(SDet.Overlap(A,Bref) == Approx(ov));
  CHECK(SDet.Overlap(Aref,Bref) == Approx(ov));

  // Test array_view
  CHECK(SDet.Overlap(A(A.extension(0),A.extension(1)),B) == Approx(ov));
  CHECK(SDet.Overlap(A,B(B.extension(0),B.extension(1))) == Approx(ov));

  array A_ = A({0,2},{0,3});
  array B_ = B({0,3},{0,2});
  REQUIRE(SDet.Overlap(A({0,2},{0,3}),
                       B({0,3},{0,2})) == Approx(ov2));
  CHECK(SDet.Overlap(A({0,2},{0,3}),B_) == Approx(ov2));
  CHECK(SDet.Overlap(A_,B({0,3},{0,2})) == Approx(ov2));


  // Density Matrices 
  vector v_ref = {
   1.3714293774508117,  -0.8858936704416944,  -0.0173026107508136,   0.1107367088052098,
   0.4324931745588929,  -0.0315365157181073,  -0.0201471975726192,   0.1289420644647630,
  -0.0999225770859535,   0.2383246559068652,   1.0046547784356801,  -0.0297905819883582,
  -1.1556702097314666,   2.7563811213471010,   0.0538355687763097,   0.6554523598316143
  };
  vector vc_ref = {
  -3.412595172195462,   7.269792211944457,   0.376363129139538,   1.591275973506948,
   0.510630761659672,   0.800731599519788,  -0.707016960946879,   0.524908550060026,
   1.072417068147838,  -3.116820833346399,   0.446937093098704,  -0.860397395831702
  };
  vector v_ref_2 = {
   0.7544916699938296,   0.6202315705419044,  -0.0193822365794349,
   0.3020580015244456,   0.2369061014119270,   0.0238466843308772,
   0.1089615621937494,  -0.2752713150157893,   1.0086022285942438
  };
  vector vc_ref_2 = {
  -0.492541105586005,  -0.860948785887990,   1.276904649559000,
   0.499074443758847,   0.581285615767123,  -0.486915175492723
  };

  array_ref g_ref(v_ref.data(),{NMO,NMO});
  array_ref gc_ref(vc_ref.data(),{NEL,NMO});
  array_ref g_ref_2(v_ref_2.data(),{3,3});
  array_ref gc_ref_2(vc_ref_2.data(),{2,3});

  array G({NMO,NMO});
  array Gc({NEL,NMO});

  ov_=SDet.MixedDensityMatrix(A,B,G,false); check(G,g_ref);
  ov_=SDet.MixedDensityMatrix(Aref,B,G,false); check(G,g_ref);
  ov_=SDet.MixedDensityMatrix(A,Bref,G,false); check(G,g_ref);
  ov_=SDet.MixedDensityMatrix(Aref,Bref,G,false); check(G,g_ref);

  ov_=SDet.MixedDensityMatrix(A({0,2},{0,3}),
                          B({0,3},{0,2}),
                          G({0,3},{0,3}),false);
  check(G({0,3},{0,3}),g_ref_2);
  ov_=SDet.MixedDensityMatrix(A_,
                          B({0,3},{0,2}),
                          G({0,3},{0,3}),false);
  check(G({0,3},{0,3}),g_ref_2);
  ov_=SDet.MixedDensityMatrix(A({0,2},{0,3}),
                          B_,
                          G({0,3},{0,3}),false);
  check(G({0,3},{0,3}),g_ref_2);

  ov_=SDet.MixedDensityMatrix(A,B,Gc,true); check(Gc,gc_ref);
  ov_=SDet.MixedDensityMatrix(Aref,B,Gc,true); check(Gc,gc_ref);
  ov_=SDet.MixedDensityMatrix(A,Bref,Gc,true); check(Gc,gc_ref);
  ov_=SDet.MixedDensityMatrix(Aref,Bref,Gc,true); check(Gc,gc_ref);

  ov_=SDet.MixedDensityMatrix(A({0,2},{0,3}),
                          B({0,3},{0,2}),
                          Gc({0,2},{0,3}),true);
  check(Gc({0,2},{0,3}),gc_ref_2);
  ov_=SDet.MixedDensityMatrix(A_,
                          B({0,3},{0,2}),
                          Gc({0,2},{0,3}),true);
  check(Gc({0,2},{0,3}),gc_ref_2);
  ov_=SDet.MixedDensityMatrix(A({0,2},{0,3}),
                          B_,
                          Gc({0,2},{0,3}),true);
  check(Gc({0,2},{0,3}),gc_ref_2);

  array Q=B;

  // Orthogonalize
  Type detR = SDet.Orthogonalize(Q);
  CHECK( ov_=SDet.Overlap_noHerm(Q,Q) == Approx(1.0)  );

}

TEST_CASE("SDetOps_double_mpi3", "[sdet_ops]")
{

  Communicate *c = OHMMS::Controller;

  using boost::mpi3::shared_communicator;
  auto world = boost::mpi3::environment::get_world_instance();
  shared_communicator node = world.split_shared(world.rank());

  const int NMO = 4;
  const int NEL = 3;

  using Type = RealType;
  using vector = std::vector<Type>;
  using array = boost::multi::array<Type,2>;
  using array_ref = boost::multi::array_ref<Type,2>;

  const Type ov = 5.10443199999999;
  const Type ov2 = -11.0204000000000;

  // some arbitrary matrices
  vector m_a = {
//   0.90000,   2.40000,   3.00000,
//   0.40000,   1.00000,   1.20000,
//   1.40000,   1.60000,   3.60000,
//   0.40000,   0.20000,   0.10000
    0.90000, 0.40000, 1.40000, 0.40000,
    2.40000, 1.00000, 1.60000, 0.20000,
    3.00000, 1.20000, 3.60000, 0.10000
  };
  vector m_b = {
   1.90000,   1.40000,   0.40000,
   1.40000,   0.20000,   2.20000,
   0.40000,   2.60000,   0.60000,
   1.10000,   0.30000,   0.90000
  };

  array A({NEL,NMO});
  array B({NMO,NEL});

  for(int i=0, k=0; i<A.size(0); i++)
    for(int j=0; j<A.size(1); j++,k++)
       A[i][j] = m_a[k];

  for(int i=0, k=0; i<B.size(0); i++)
    for(int j=0; j<B.size(1); j++,k++)
       B[i][j] = m_b[k];

  array_ref Aref(m_a.data(),{NEL,NMO});
  array_ref Bref(m_b.data(),{NMO,NEL});

  SlaterDetOperations SDet( SlaterDetOperations_shared<Type>(NMO,NEL) );

  // Overlaps 
  CHECK(SDet.Overlap(A,B,node) == Approx(ov));
  CHECK(SDet.Overlap(Aref,B,node) == Approx(ov));
  CHECK(SDet.Overlap(A,Bref,node) == Approx(ov));
  CHECK(SDet.Overlap(Aref,Bref,node) == Approx(ov));

  // Test array_view
  CHECK(SDet.Overlap(A(A.extension(0),A.extension(1)),B,node) == Approx(ov));
  CHECK(SDet.Overlap(A,B(B.extension(0),B.extension(1)),node) == Approx(ov));

  array A_ = A({0,2},{0,3});
  array B_ = B({0,3},{0,2});
  REQUIRE(SDet.Overlap(A({0,2},{0,3}),
                       B({0,3},{0,2}),node) == Approx(ov2));
  CHECK(SDet.Overlap(A({0,2},{0,3}),B_) == Approx(ov2));
  CHECK(SDet.Overlap(A_,B({0,3},{0,2}),node) == Approx(ov2));

  shared_communicator node_ = node.split(node.rank()%2,node.rank());
  CHECK(SDet.Overlap(A,B,node_) == Approx(ov));
  REQUIRE(SDet.Overlap(A({0,2},{0,3}),
                       B({0,3},{0,2}),node_) == Approx(ov2));


  // Density Matrices
  vector v_ref = {
   1.3714293774508117,  -0.8858936704416944,  -0.0173026107508136,   0.1107367088052098,
   0.4324931745588929,  -0.0315365157181073,  -0.0201471975726192,   0.1289420644647630,
  -0.0999225770859535,   0.2383246559068652,   1.0046547784356801,  -0.0297905819883582,
  -1.1556702097314666,   2.7563811213471010,   0.0538355687763097,   0.6554523598316143
  };
  vector vc_ref = {
  -3.412595172195462,   7.269792211944457,   0.376363129139538,   1.591275973506948,
   0.510630761659672,   0.800731599519788,  -0.707016960946879,   0.524908550060026,
   1.072417068147838,  -3.116820833346399,   0.446937093098704,  -0.860397395831702
  };
  vector v_ref_2 = {
   0.7544916699938296,   0.6202315705419044,  -0.0193822365794349,
   0.3020580015244456,   0.2369061014119270,   0.0238466843308772,
   0.1089615621937494,  -0.2752713150157893,   1.0086022285942438
  };
  vector vc_ref_2 = {
  -0.492541105586005,  -0.860948785887990,   1.276904649559000,
   0.499074443758847,   0.581285615767123,  -0.486915175492723
  };

  array_ref g_ref(v_ref.data(),{NMO,NMO});
  array_ref gc_ref(vc_ref.data(),{NEL,NMO});
  array_ref g_ref_2(v_ref_2.data(),{3,3});
  array_ref gc_ref_2(vc_ref_2.data(),{2,3});

  boost::multi::array<Type,1,shared_allocator<Type>> SMbuff(iextensions<1u>{NMO*(NMO+NEL)},
                                                            shared_allocator<Type>{node});  

  array_ref G(to_address(SMbuff.origin()),{NMO,NMO});
  array_ref Gc(to_address(SMbuff.origin())+NMO*NMO,{NEL,NMO});

  ov_=SDet.MixedDensityMatrix(A,B,G,node,false); check(G,g_ref);
  ov_=SDet.MixedDensityMatrix(Aref,B,G,node,false); check(G,g_ref);
  ov_=SDet.MixedDensityMatrix(A,Bref,G,node,false); check(G,g_ref);
  ov_=SDet.MixedDensityMatrix(Aref,Bref,G,node,false); check(G,g_ref);

  ov_=SDet.MixedDensityMatrix(A({0,2},{0,3}),
                          B({0,3},{0,2}),
                          G({0,3},{0,3}),node,false);
  check(G({0,3},{0,3}),g_ref_2);
  ov_=SDet.MixedDensityMatrix(A_,
                          B({0,3},{0,2}),
                          G({0,3},{0,3}),node,false);
  check(G({0,3},{0,3}),g_ref_2);
  ov_=SDet.MixedDensityMatrix(A({0,2},{0,3}),
                          B_,
                          G({0,3},{0,3}),node,false);
  check(G({0,3},{0,3}),g_ref_2);

  ov_=SDet.MixedDensityMatrix(A,B,Gc,node,true); check(Gc,gc_ref);
  ov_=SDet.MixedDensityMatrix(Aref,B,Gc,node,true); check(Gc,gc_ref);
  ov_=SDet.MixedDensityMatrix(A,Bref,Gc,node,true); check(Gc,gc_ref);
  ov_=SDet.MixedDensityMatrix(Aref,Bref,Gc,node,true); check(Gc,gc_ref);

  ov_=SDet.MixedDensityMatrix(A({0,2},{0,3}),
                          B({0,3},{0,2}),
                          Gc({0,2},{0,3}),node,true);
  check(Gc({0,2},{0,3}),gc_ref_2);
  ov_=SDet.MixedDensityMatrix(A_,
                          B({0,3},{0,2}),
                          Gc({0,2},{0,3}),node,true);
  check(Gc({0,2},{0,3}),gc_ref_2);
  ov_=SDet.MixedDensityMatrix(A({0,2},{0,3}),
                          B_,
                          Gc({0,2},{0,3}),node,true);
  check(Gc({0,2},{0,3}),gc_ref_2);

  boost::multi::array<Type,1,shared_allocator<Type>> SMbuff2(iextensions<1u>{NMO*(NMO+NEL)},
                                                            shared_allocator<Type>{node_});

  array_ref G2(to_address(SMbuff2.origin()),{NMO,NMO});
  array_ref Gc2(to_address(SMbuff2.origin())+NMO*NMO,{NEL,NMO});

  // switch comm
  ov_=SDet.MixedDensityMatrix(A,B,G2,node_,false); check(G2,g_ref);
  ov_=SDet.MixedDensityMatrix(A({0,2},{0,3}),
                          B({0,3},{0,2}),
                          G2({0,3},{0,3}),node_,false);
  check(G2({0,3},{0,3}),g_ref_2);
  ov_=SDet.MixedDensityMatrix(A,B,Gc2,node_,true); check(Gc2,gc_ref);
  ov_=SDet.MixedDensityMatrix(A({0,2},{0,3}),
                          B({0,3},{0,2}),
                          Gc2({0,2},{0,3}),node_,true);
  check(Gc2({0,2},{0,3}),gc_ref_2);

}
*/

template<class Allocator, class BufferManager>
void SDetOps_complex_serial(Allocator alloc, BufferManager b)
{
  static_assert(std::is_same<typename Allocator::value_type, ComplexType>::value, "Incorrect type.\n");

  const int NMO = 4;
  const int NEL = 3;

  using Type      = ComplexType;
  using vector    = std::vector<Type>;
  using array     = boost::multi::array<Type, 2, Allocator>;
  using array_ref = boost::multi::array_ref<Type, 2, typename Allocator::pointer>;
  using array_ptr = boost::multi::array_ptr<Type, 2, typename Allocator::pointer>;
  using namespace std::complex_literals;

  const Type ov  = -7.62332599999999 + 22.20453200000000i;
  const Type ov2 = -10.37150000000000 - 7.15750000000000i;

  // some arbitrary matrices
  vector m_a = {0.90000 + 0.10000i, 0.40000 + 0.40000i, 1.40000 + 0.20000i, 0.40000 + 0.50000i,
                2.40000 + 0.20000i, 1.00000 + 0.50000i, 1.60000 + 0.30000i, 0.20000 + 0.10000i,
                3.00000 + 0.30000i, 1.20000 + 0.10000i, 3.60000 + 0.40000i, 0.10000 + 0.20000i};
  vector m_b = {1.90000 + 0.60000i, 1.40000 + 0.70000i, 0.40000 + 0.80000i, 1.40000 + 0.90000i,
                0.20000 + 0.50000i, 2.20000 + 0.60000i, 0.40000 + 0.70000i, 2.60000 + 0.80000i,
                0.60000 + 0.90000i, 1.10000 + 0.50000i, 0.30000 + 0.60000i, 0.90000 + 0.70000i};

  array A({NEL, NMO}, alloc);
  copy_n(m_a.data(), m_a.size(), A.origin());
  array B({NMO, NEL}, alloc);
  copy_n(m_b.data(), m_b.size(), B.origin());

  array_ref Aref(A.origin(), {NEL, NMO});
  array_ref Bref(B.origin(), {NMO, NEL});

  SlaterDetOperations SDet(SlaterDetOperations_serial<Type, BufferManager>(NMO, NEL, b));

  /**** Overlaps ****/
  //SECTION("Overlaps")
  {
    Type ov_;
    ov_ = SDet.Overlap(A, B, 0.0);
    myCHECK(ov_, ov);
    ov_ = SDet.Overlap(Aref, B, 0.0);
    myCHECK(ov_, ov);
    ov_ = SDet.Overlap(A, Bref, 0.0);
    myCHECK(ov_, ov);
    ov_ = SDet.Overlap(Aref, Bref, 0.0);
    myCHECK(ov_, ov);
  }

  // Test array_view
  //SECTION("array_view")
  {
    Type ov_;
    ov_ = SDet.Overlap(A(A.extension(0), A.extension(1)), B, 0.0);
    myCHECK(ov_, ov);
    ov_ = SDet.Overlap(A, B(B.extension(0), B.extension(1)), 0.0);
    myCHECK(ov_, ov);
  }

  // copy not yet working with device_pointer
  //SECTION("copy_overlap")
  {
    array A_ = A({0, 2}, {0, 3});
    array B_ = B({0, 3}, {0, 2});
    Type ov_;
    ov_ = SDet.Overlap(A({0, 2}, {0, 3}), B({0, 3}, {0, 2}), 0.0);
    myCHECK(ov_, ov2);
    ov_ = SDet.Overlap(A({0, 2}, {0, 3}), B_, 0.0);
    myCHECK(ov_, ov2);
    ov_ = SDet.Overlap(A_, B({0, 3}, {0, 2}), 0.0);
    myCHECK(ov_, ov2);
  }

  /**** Density Matrices *****/
  vector v_ref    = {1.17573619385025996 - 0.01580426445014660i,  -0.25295981756593167 + 0.28594469607401085i,
                  -0.07724502823533341 - 0.09687959052155870i, 0.30512858581808422 - 0.04506898328729603i,

                  0.17912592806889663 + 0.08374315906672802i,  0.59381451118048767 + 0.13438888951771200i,
                  -0.02021475320201610 - 0.13737561982083193i, 0.32095003919745313 + 0.12832750636154097i,

                  -0.04919549564646425 + 0.05402065222741825i, -0.00286990878355775 - 0.15806420733175885i,
                  1.05069081188635494 + 0.00793429988912356i,  -0.08048239150997794 + 0.09917405634760490i,

                  -0.46434598219794548 - 0.09896890422706731i, 0.87746748807670427 - 0.53417787485950319i,
                  0.12162735438647077 + 0.31042401735800573i,  0.17975848308289613 - 0.12651892495668898i};
  vector vc_ref   = {-0.8721879971495297 + 0.6787593377585239i, 0.3244278932250768 - 1.8083537898275881i,
                   0.6130713546530860 + 0.0399736955598931i,  0.0132562806444336 - 0.2882495766584950i,

                   0.8853626557603183 + 0.0978868569224204i,  -0.0598704127345155 - 0.0470889603064014i,
                   -0.7392693424168300 - 0.0715317395994149i, 0.3721269544963505 - 0.1797896522886788i,

                   -0.0567190307022984 - 0.3114847157576828i, -0.1290126440468128 + 0.6815705660308808i,
                   0.3787012855335005 + 0.0039188686237135i,  -0.2005543456941538 + 0.2100886953142371i};
  vector v_ref_2  = {0.7361983496013835 - 0.0956505507662245i, 0.6467449689807925 + 0.2297471806893873i,
                    0.0189270005620390 - 0.1727975708935829i,

                    0.2986893588843604 - 0.0099955730815030i, 0.2696068479051557 + 0.0292039710860386i,
                    0.0497734835066391 + 0.1783200397050796i,

                    0.1020030246826092 + 0.0344707468383766i, -0.2500340021988402 - 0.0826863644855427i,
                    0.9941948024934623 + 0.0664465796801866i};
  vector vc_ref_2 = {-0.489369975192701 + 0.103038673040713i, -0.858850485405126 - 0.275734238124941i,
                     1.219791948842170 + 0.118626922447301i,  0.486337033653898 - 0.098631569047656i,
                     0.595497821653010 + 0.185288949671560i,  -0.455373025165330 - 0.129360996228044i};

  boost::multi::array_ref<Type, 2> g_ref(v_ref.data(), {NMO, NMO});
  boost::multi::array_ref<Type, 2> gc_ref(vc_ref.data(), {NEL, NMO});
  boost::multi::array_ref<Type, 2> g_ref_2(v_ref_2.data(), {3, 3});
  boost::multi::array_ref<Type, 2> gc_ref_2(vc_ref_2.data(), {2, 3});

  array G({NMO, NMO}, alloc);
  array Gc({NEL, NMO}, alloc);
  //SECTION("density_matrices")
  {
    Type ov_;
    array A_ = A({0, 2}, {0, 3});
    array B_ = B({0, 3}, {0, 2});
    ov_      = SDet.MixedDensityMatrix(A, B, G, 0.0, false);
    check(G, g_ref);
    ov_ = SDet.MixedDensityMatrix(Aref, B, G, 0.0, false);
    check(G, g_ref);
    ov_ = SDet.MixedDensityMatrix(A, Bref, G, 0.0, false);
    check(G, g_ref);
    ov_ = SDet.MixedDensityMatrix(Aref, Bref, G, 0.0, false);
    check(G, g_ref);

    ov_ = SDet.MixedDensityMatrix(A({0, 2}, {0, 3}), B({0, 3}, {0, 2}), G({0, 3}, {0, 3}), 0.0, false);
    check(G({0, 3}, {0, 3}), g_ref_2);
    ov_ = SDet.MixedDensityMatrix(A_, B({0, 3}, {0, 2}), G({0, 3}, {0, 3}), 0.0, false);
    check(G({0, 3}, {0, 3}), g_ref_2);
    ov_ = SDet.MixedDensityMatrix(A({0, 2}, {0, 3}), B_, G({0, 3}, {0, 3}), 0.0, false);
    check(G({0, 3}, {0, 3}), g_ref_2);

    ov_ = SDet.MixedDensityMatrix(A, B, Gc, 0.0, true);
    check(Gc, gc_ref);
    ov_ = SDet.MixedDensityMatrix(Aref, B, Gc, 0.0, true);
    check(Gc, gc_ref);
    ov_ = SDet.MixedDensityMatrix(A, Bref, Gc, 0.0, true);
    check(Gc, gc_ref);
    ov_ = SDet.MixedDensityMatrix(Aref, Bref, Gc, 0.0, true);
    check(Gc, gc_ref);
    myCHECK(ov_, ov);

    ov_ = SDet.MixedDensityMatrix(A({0, 2}, {0, 3}), B({0, 3}, {0, 2}), Gc({0, 2}, {0, 3}), 0.0, true);
    check(Gc({0, 2}, {0, 3}), gc_ref_2);
    ov_ = SDet.MixedDensityMatrix(A_, B({0, 3}, {0, 2}), Gc({0, 2}, {0, 3}), 0.0, true);
    check(Gc({0, 2}, {0, 3}), gc_ref_2);
    ov_ = SDet.MixedDensityMatrix(A({0, 2}, {0, 3}), B_, Gc({0, 2}, {0, 3}), 0.0, true);
    check(Gc({0, 2}, {0, 3}), gc_ref_2);
  }
  // Orthogonalize
  //SECTION("orthogonalize")
  {
    array Q = B;
    SDet.Orthogonalize(Q, 0.0);
    Type ov_ = SDet.Overlap_noHerm(Q, Q, 0.0);
    myCHECK(ov_, std::complex<double>(1., 0.));
  }

  // Batched
  // TODO fix CPU.
#if defined(ENABLE_CUDA) || defined(ENABLE_HIP)
  //SECTION("batched_density_matrix")
  {
    boost::multi::array<Type, 3, Allocator> Gw({3, NMO, NMO}, alloc);
    //std::vector<array_ptr> RA;
    std::vector<decltype(&Aref)> RA;
    std::vector<decltype(&Bref)> RB;
    std::vector<decltype(&Gw[0])> Gwv;
    RA.reserve(3);
    RB.reserve(3);
    Gwv.reserve(3);
    boost::multi::array<Type, 1, Allocator> ovlp(iextensions<1u>{3});
    for (int i = 0; i < 3; i++)
    {
      RA.emplace_back(&Aref);
      RB.emplace_back(&Bref);
      Gwv.emplace_back(&Gw[i]);
    }
    Type log_ovlp;
    Type ov_ref = -7.623325999999989 + 22.20453200000001i;

    SDet.BatchedDensityMatrices(RA, RB, Gwv, log_ovlp, ovlp, false);
    //SECTION("greens_function")
    {
      for (int i = 0; i < 3; i++)
      {
        check(*Gwv[i], g_ref);
      }
    }
    //SECTION("overlap")
    {
      for (int i = 0; i < 3; i++)
      {
        myCHECK(ovlp[i], ov_ref);
      }
    }
    SDet.BatchedMixedDensityMatrix(RA, RB, Gw, log_ovlp, ovlp, false);
    {
      for (int i = 0; i < 3; i++)
      {
        myCHECK(ovlp[i], ov_ref);
      }
    }
    {
      for (int i = 0; i < 3; i++)
      {
        check(Gw[i], g_ref);
      }
    }
  }
  //SECTION("batched_mixed_density_matrix")
  {
    boost::multi::array<Type, 3, Allocator> Gw({3, NEL, NMO}, alloc);
    //std::vector<array_ptr> RA;
    std::vector<decltype(&Aref)> RA;
    std::vector<decltype(&Bref)> RB;
    RA.reserve(3);
    RB.reserve(3);
    boost::multi::array<Type, 1, Allocator> ovlp(iextensions<1u>{3});
    for (int i = 0; i < 3; i++)
    {
      RA.emplace_back(&Aref);
      RB.emplace_back(&Bref);
    }
    Type log_ovlp;
    Type ov_ref = -7.623325999999989 + 22.20453200000001i;

    SDet.BatchedMixedDensityMatrix(RA, RB, Gw, log_ovlp, ovlp, true);
    //SECTION("greens_function")
    //SECTION("overlap")
    {
      for (int i = 0; i < 3; i++)
      {
        myCHECK(ovlp[i], ov_ref);
      }
    }
    {
      for (int i = 0; i < 3; i++)
      {
        check(Gw[i], gc_ref);
      }
    }
  }
#endif
}

TEST_CASE("SDetOps_complex_mpi3", "[sdet_ops]")
{
  Communicate* c = OHMMS::Controller;

  using boost::mpi3::shared_communicator;
  auto world               = boost::mpi3::environment::get_world_instance();
  shared_communicator node = world.split_shared(world.rank());
  setup_memory_managers(node, 10uL * 1024uL * 1024uL);

  const int NMO = 4;
  const int NEL = 3;

  using Type      = ComplexType;
  using vector    = std::vector<Type>;
  using array     = boost::multi::array<Type, 2>;
  using array_ref = boost::multi::array_ref<Type, 2>;
  using namespace std::complex_literals;

  const Type ov  = -7.62332599999999 + 22.20453200000000i;
  const Type ov2 = -10.37150000000000 - 7.15750000000000i;

  // some arbitrary matrices
  vector m_a = {//   0.90000 + 0.10000i,   2.40000 + 0.20000i,   3.00000 + 0.30000i,
                //   0.40000 + 0.40000i,   1.00000 + 0.50000i,   1.20000 + 0.10000i,
                //   1.40000 + 0.20000i,   1.60000 + 0.30000i,   3.60000 + 0.40000i,
                //   0.40000 + 0.50000i,   0.20000 + 0.10000i,   0.10000 + 0.20000i
                0.90000 + 0.10000i, 0.40000 + 0.40000i, 1.40000 + 0.20000i, 0.40000 + 0.50000i,
                2.40000 + 0.20000i, 1.00000 + 0.50000i, 1.60000 + 0.30000i, 0.20000 + 0.10000i,
                3.00000 + 0.30000i, 1.20000 + 0.10000i, 3.60000 + 0.40000i, 0.10000 + 0.20000i};
  vector m_b = {1.90000 + 0.60000i, 1.40000 + 0.70000i, 0.40000 + 0.80000i, 1.40000 + 0.90000i,
                0.20000 + 0.50000i, 2.20000 + 0.60000i, 0.40000 + 0.70000i, 2.60000 + 0.80000i,
                0.60000 + 0.90000i, 1.10000 + 0.50000i, 0.30000 + 0.60000i, 0.90000 + 0.70000i};

  array A({NEL, NMO});
  array B({NMO, NEL});

  for (int i = 0, k = 0; i < std::get<0>(A.sizes()); i++)
    for (int j = 0; j < std::get<1>(A.sizes()); j++, k++)
      A[i][j] = m_a[k];

  for (int i = 0, k = 0; i < std::get<0>(B.sizes()); i++)
    for (int j = 0; j < std::get<1>(B.sizes()); j++, k++)
      B[i][j] = m_b[k];

  array_ref Aref(m_a.data(), {NEL, NMO});
  array_ref Bref(m_b.data(), {NMO, NEL});

  //SlaterDetOperations SDet( SlaterDetOperations_shared<Type>(NMO,NEL) );
  SlaterDetOperations_shared<ComplexType> SDet(NMO, NEL);

  /**** Overlaps ****/
  Type ov_;
  ov_ = SDet.Overlap(A, B, 0.0, node);
  myCHECK(ov_, ov);
  ov_ = SDet.Overlap(Aref, B, 0.0, node);
  myCHECK(ov_, ov);
  ov_ = SDet.Overlap(A, Bref, 0.0, node);
  myCHECK(ov_, ov);
  ov_ = SDet.Overlap(Aref, Bref, 0.0, node);
  myCHECK(ov_, ov);

  // Test array_view
  ov_ = SDet.Overlap(A(A.extension(0), A.extension(1)), B, 0.0, node);
  myCHECK(ov_, ov);
  ov_ = SDet.Overlap(A, B(B.extension(0), B.extension(1)), 0.0, node);
  myCHECK(ov_, ov);

  array A_ = A({0, 2}, {0, 3});
  array B_ = B({0, 3}, {0, 2});
  ov_      = SDet.Overlap(A({0, 2}, {0, 3}), B({0, 3}, {0, 2}), 0.0, node);
  myCHECK(ov_, ov2);
  ov_ = SDet.Overlap(A({0, 2}, {0, 3}), B_, 0.0);
  myCHECK(ov_, ov2);
  ov_ = SDet.Overlap(A_, B({0, 3}, {0, 2}), 0.0, node);
  myCHECK(ov_, ov2);

  shared_communicator node_ = node.split(node.rank() % 2, node.rank());
  ov_                       = SDet.Overlap(A, B, 0.0, node_);
  myCHECK(ov_, ov);
  ov_ = SDet.Overlap(A({0, 2}, {0, 3}), B({0, 3}, {0, 2}), 0.0, node_);
  myCHECK(ov_, ov2);

  /**** Density Matrices *****/
  vector v_ref    = {1.17573619385025996 - 0.01580426445014660i,  -0.25295981756593167 + 0.28594469607401085i,
                  -0.07724502823533341 - 0.09687959052155870i, 0.30512858581808422 - 0.04506898328729603i,

                  0.17912592806889663 + 0.08374315906672802i,  0.59381451118048767 + 0.13438888951771200i,
                  -0.02021475320201610 - 0.13737561982083193i, 0.32095003919745313 + 0.12832750636154097i,

                  -0.04919549564646425 + 0.05402065222741825i, -0.00286990878355775 - 0.15806420733175885i,
                  1.05069081188635494 + 0.00793429988912356i,  -0.08048239150997794 + 0.09917405634760490i,

                  -0.46434598219794548 - 0.09896890422706731i, 0.87746748807670427 - 0.53417787485950319i,
                  0.12162735438647077 + 0.31042401735800573i,  0.17975848308289613 - 0.12651892495668898i};
  vector vc_ref   = {-0.8721879971495297 + 0.6787593377585239i, 0.3244278932250768 - 1.8083537898275881i,
                   0.6130713546530860 + 0.0399736955598931i,  0.0132562806444336 - 0.2882495766584950i,

                   0.8853626557603183 + 0.0978868569224204i,  -0.0598704127345155 - 0.0470889603064014i,
                   -0.7392693424168300 - 0.0715317395994149i, 0.3721269544963505 - 0.1797896522886788i,

                   -0.0567190307022984 - 0.3114847157576828i, -0.1290126440468128 + 0.6815705660308808i,
                   0.3787012855335005 + 0.0039188686237135i,  -0.2005543456941538 + 0.2100886953142371i};
  vector v_ref_2  = {0.7361983496013835 - 0.0956505507662245i, 0.6467449689807925 + 0.2297471806893873i,
                    0.0189270005620390 - 0.1727975708935829i,

                    0.2986893588843604 - 0.0099955730815030i, 0.2696068479051557 + 0.0292039710860386i,
                    0.0497734835066391 + 0.1783200397050796i,

                    0.1020030246826092 + 0.0344707468383766i, -0.2500340021988402 - 0.0826863644855427i,
                    0.9941948024934623 + 0.0664465796801866i};
  vector vc_ref_2 = {-0.489369975192701 + 0.103038673040713i, -0.858850485405126 - 0.275734238124941i,
                     1.219791948842170 + 0.118626922447301i,  0.486337033653898 - 0.098631569047656i,
                     0.595497821653010 + 0.185288949671560i,  -0.455373025165330 - 0.129360996228044i};

  boost::multi::array_ref<Type, 2> g_ref(v_ref.data(), {NMO, NMO});
  boost::multi::array_ref<Type, 2> gc_ref(vc_ref.data(), {NEL, NMO});
  boost::multi::array_ref<Type, 2> g_ref_2(v_ref_2.data(), {3, 3});
  boost::multi::array_ref<Type, 2> gc_ref_2(vc_ref_2.data(), {2, 3});

  boost::multi::array<Type, 1, shared_allocator<Type>> SMbuff(iextensions<1u>{NMO * (NMO + NEL)},
                                                              shared_allocator<Type>{node});

  array_ref G(to_address(SMbuff.origin()), {NMO, NMO});
  array_ref Gc(to_address(SMbuff.origin()) + NMO * NMO, {NEL, NMO});

  ov_ = SDet.MixedDensityMatrix(A, B, G, 0.0, node, false);
  check(G, g_ref);
  ov_ = SDet.MixedDensityMatrix(Aref, B, G, 0.0, node, false);
  check(G, g_ref);
  ov_ = SDet.MixedDensityMatrix(A, Bref, G, 0.0, node, false);
  check(G, g_ref);
  ov_ = SDet.MixedDensityMatrix(Aref, Bref, G, 0.0, node, false);
  check(G, g_ref);

  ov_ = SDet.MixedDensityMatrix(A({0, 2}, {0, 3}), B({0, 3}, {0, 2}), G({0, 3}, {0, 3}), 0.0, node, false);
  check(G({0, 3}, {0, 3}), g_ref_2);
  ov_ = SDet.MixedDensityMatrix(A_, B({0, 3}, {0, 2}), G({0, 3}, {0, 3}), 0.0, node, false);
  check(G({0, 3}, {0, 3}), g_ref_2);
  ov_ = SDet.MixedDensityMatrix(A({0, 2}, {0, 3}), B_, G({0, 3}, {0, 3}), 0.0, node, false);
  check(G({0, 3}, {0, 3}), g_ref_2);

  ov_ = SDet.MixedDensityMatrix(A, B, Gc, 0.0, node, true);
  check(Gc, gc_ref);
  ov_ = SDet.MixedDensityMatrix(Aref, B, Gc, 0.0, node, true);
  check(Gc, gc_ref);
  ov_ = SDet.MixedDensityMatrix(A, Bref, Gc, 0.0, node, true);
  check(Gc, gc_ref);
  ov_ = SDet.MixedDensityMatrix(Aref, Bref, Gc, 0.0, node, true);
  check(Gc, gc_ref);

  ov_ = SDet.MixedDensityMatrix(A({0, 2}, {0, 3}), B({0, 3}, {0, 2}), Gc({0, 2}, {0, 3}), 0.0, node, true);
  check(Gc({0, 2}, {0, 3}), gc_ref_2);
  ov_ = SDet.MixedDensityMatrix(A_, B({0, 3}, {0, 2}), Gc({0, 2}, {0, 3}), 0.0, node, true);
  check(Gc({0, 2}, {0, 3}), gc_ref_2);
  ov_ = SDet.MixedDensityMatrix(A({0, 2}, {0, 3}), B_, Gc({0, 2}, {0, 3}), 0.0, node, true);
  check(Gc({0, 2}, {0, 3}), gc_ref_2);

  boost::multi::array<Type, 1, shared_allocator<Type>> SMbuff2(iextensions<1u>{NMO * (NMO + NEL)},
                                                               shared_allocator<Type>{node_});

  array_ref G2(to_address(SMbuff2.origin()), {NMO, NMO});
  array_ref Gc2(to_address(SMbuff2.origin()) + NMO * NMO, {NEL, NMO});

  // switch comm
  ov_ = SDet.MixedDensityMatrix(A, B, G2, 0.0, node_, false);
  check(G2, g_ref);
  ov_ = SDet.MixedDensityMatrix(A({0, 2}, {0, 3}), B({0, 3}, {0, 2}), G2({0, 3}, {0, 3}), 0.0, node_, false);
  check(G2({0, 3}, {0, 3}), g_ref_2);
  ov_ = SDet.MixedDensityMatrix(A, B, Gc2, 0.0, node_, true);
  check(Gc2, gc_ref);
  ov_ = SDet.MixedDensityMatrix(A({0, 2}, {0, 3}), B({0, 3}, {0, 2}), Gc2({0, 2}, {0, 3}), 0.0, node_, true);
  check(Gc2({0, 2}, {0, 3}), gc_ref_2);
  release_memory_managers();
}

TEST_CASE("SDetOps_complex_csr", "[sdet_ops]")
{
  Communicate* c = OHMMS::Controller;

  using boost::mpi3::shared_communicator;

  auto world               = boost::mpi3::environment::get_world_instance();
  shared_communicator node = world.split_shared(world.rank());
  setup_memory_managers(node, 10uL * 1024uL * 1024uL);

  const int NMO = 4;
  const int NEL = 3;

  using Type      = std::complex<double>;
  using vector    = std::vector<Type>;
  using array     = boost::multi::array<Type, 2>;
  using array_ref = boost::multi::array_ref<Type, 2>;
  using namespace std::complex_literals;
  using csr_matrix = ma::sparse::csr_matrix<Type, int, int, shared_allocator<Type>, ma::sparse::is_root>;

  const Type ov  = -7.62332599999999 + 22.20453200000000i;
  const Type ov2 = -10.37150000000000 - 7.15750000000000i;

  // some arbitrary matrices
  vector m_a = {0.90000 + 0.10000i, 2.40000 + 0.20000i, 3.00000 + 0.30000i, 0.40000 + 0.40000i,
                1.00000 + 0.50000i, 1.20000 + 0.10000i, 1.40000 + 0.20000i, 1.60000 + 0.30000i,
                3.60000 + 0.40000i, 0.40000 + 0.50000i, 0.20000 + 0.10000i, 0.10000 + 0.20000i};
  vector m_b = {1.90000 + 0.60000i, 1.40000 + 0.70000i, 0.40000 + 0.80000i, 1.40000 + 0.90000i,
                0.20000 + 0.50000i, 2.20000 + 0.60000i, 0.40000 + 0.70000i, 2.60000 + 0.80000i,
                0.60000 + 0.90000i, 1.10000 + 0.50000i, 0.30000 + 0.60000i, 0.90000 + 0.70000i};

  array A({NMO, NEL}); // Will be transposed when Acsr is built
  array B({NMO, NEL});

  for (int i = 0, k = 0; i < std::get<0>(A.sizes()); i++)
    for (int j = 0; j < std::get<1>(A.sizes()); j++, k++)
      A[i][j] = m_a[k];

  for (int i = 0, k = 0; i < std::get<0>(B.sizes()); i++)
    for (int j = 0; j < std::get<1>(B.sizes()); j++, k++)
      B[i][j] = m_b[k];

  boost::multi::array_ref<Type, 2> Bref(m_b.data(), {NMO, NEL});

  csr_matrix Acsr(csr::shm::construct_csr_matrix_single_input<csr_matrix>(A, 0.0, 'T', node));

  //SlaterDetOperations SDet( SlaterDetOperations_shared<Type>(NMO,NEL) );
  SlaterDetOperations_shared<ComplexType> SDet(NMO, NEL);

  /**** Overlaps ****/
  Type ov_;
  ov_ = SDet.Overlap(Acsr, B, 0.0, node);
  myCHECK(ov_, ov);
  ov_ = SDet.Overlap(Acsr, Bref, 0.0, node);
  myCHECK(ov_, ov);

  ov_ = SDet.Overlap(Acsr, B, 0.0);
  myCHECK(ov_, ov);
  ov_ = SDet.Overlap(Acsr, Bref, 0.0);
  myCHECK(ov_, ov);

  // Test array_view
  ov_ = SDet.Overlap(Acsr, B(B.extension(0), B.extension(1)), 0.0, node);
  myCHECK(ov_, ov);
  ov_ = SDet.Overlap(Acsr, B(B.extension(0), B.extension(1)), 0.0);
  myCHECK(ov_, ov);

  shared_communicator node_ = node.split(node.rank() % 2, node.rank());
  ov_                       = SDet.Overlap(Acsr, B, 0.0, node_);
  myCHECK(ov_, ov);

  array B_ = B({0, 3}, {0, 2});

  ov_ = SDet.Overlap(Acsr[{0, 2, 0, 3}], B_, 0.0);
  myCHECK(ov_, ov2);

  /**** Density Matrices *****/
  vector v_ref    = {1.17573619385025996 - 0.01580426445014660i,  -0.25295981756593167 + 0.28594469607401085i,
                  -0.07724502823533341 - 0.09687959052155870i, 0.30512858581808422 - 0.04506898328729603i,

                  0.17912592806889663 + 0.08374315906672802i,  0.59381451118048767 + 0.13438888951771200i,
                  -0.02021475320201610 - 0.13737561982083193i, 0.32095003919745313 + 0.12832750636154097i,

                  -0.04919549564646425 + 0.05402065222741825i, -0.00286990878355775 - 0.15806420733175885i,
                  1.05069081188635494 + 0.00793429988912356i,  -0.08048239150997794 + 0.09917405634760490i,

                  -0.46434598219794548 - 0.09896890422706731i, 0.87746748807670427 - 0.53417787485950319i,
                  0.12162735438647077 + 0.31042401735800573i,  0.17975848308289613 - 0.12651892495668898i};
  vector vc_ref   = {-0.8721879971495297 + 0.6787593377585239i, 0.3244278932250768 - 1.8083537898275881i,
                   0.6130713546530860 + 0.0399736955598931i,  0.0132562806444336 - 0.2882495766584950i,

                   0.8853626557603183 + 0.0978868569224204i,  -0.0598704127345155 - 0.0470889603064014i,
                   -0.7392693424168300 - 0.0715317395994149i, 0.3721269544963505 - 0.1797896522886788i,

                   -0.0567190307022984 - 0.3114847157576828i, -0.1290126440468128 + 0.6815705660308808i,
                   0.3787012855335005 + 0.0039188686237135i,  -0.2005543456941538 + 0.2100886953142371i};
  vector v_ref_2  = {0.7361983496013835 - 0.0956505507662245i, 0.6467449689807925 + 0.2297471806893873i,
                    0.0189270005620390 - 0.1727975708935829i,

                    0.2986893588843604 - 0.0099955730815030i, 0.2696068479051557 + 0.0292039710860386i,
                    0.0497734835066391 + 0.1783200397050796i,

                    0.1020030246826092 + 0.0344707468383766i, -0.2500340021988402 - 0.0826863644855427i,
                    0.9941948024934623 + 0.0664465796801866i};
  vector vc_ref_2 = {-0.489369975192701 + 0.103038673040713i, -0.858850485405126 - 0.275734238124941i,
                     1.219791948842170 + 0.118626922447301i,  0.486337033653898 - 0.098631569047656i,
                     0.595497821653010 + 0.185288949671560i,  -0.455373025165330 - 0.129360996228044i};

  boost::multi::array_ref<Type, 2> g_ref(v_ref.data(), {NMO, NMO});
  boost::multi::array_ref<Type, 2> gc_ref(vc_ref.data(), {NEL, NMO});
  boost::multi::array_ref<Type, 2> g_ref_2(v_ref_2.data(), {3, 3});
  boost::multi::array_ref<Type, 2> gc_ref_2(vc_ref_2.data(), {2, 3});

  boost::multi::array<Type, 1, shared_allocator<Type>> SMbuff(iextensions<1u>{NMO * (NMO + NEL)},
                                                              shared_allocator<Type>{node});

  array_ref G(to_address(SMbuff.origin()), {NMO, NMO});
  array_ref Gc(to_address(SMbuff.origin()) + NMO * NMO, {NEL, NMO});

  ov_ = SDet.MixedDensityMatrix(Acsr, B, G, 0.0, node, false);
  check(G, g_ref);
  ov_ = SDet.MixedDensityMatrix(Acsr, Bref, G, 0.0, node, false);
  check(G, g_ref);

  ov_ = SDet.MixedDensityMatrix(Acsr[{0, 2, 0, 3}], B({0, 3}, {0, 2}), G({0, 3}, {0, 3}), 0.0, node, false);
  check(G({0, 3}, {0, 3}), g_ref_2);
  ov_ = SDet.MixedDensityMatrix(Acsr[{0, 2, 0, 3}], B_, G({0, 3}, {0, 3}), 0.0, node, false);
  check(G({0, 3}, {0, 3}), g_ref_2);

  ov_ = SDet.MixedDensityMatrix(Acsr, B, Gc, 0.0, node, true);
  check(Gc, gc_ref);
  ov_ = SDet.MixedDensityMatrix(Acsr, Bref, Gc, 0.0, node, true);
  check(Gc, gc_ref);

  ov_ = SDet.MixedDensityMatrix(Acsr[{0, 2, 0, 3}], B({0, 3}, {0, 2}), Gc({0, 2}, {0, 3}), 0.0, node, true);
  check(Gc({0, 2}, {0, 3}), gc_ref_2);
  ov_ = SDet.MixedDensityMatrix(Acsr[{0, 2, 0, 3}], B_, Gc({0, 2}, {0, 3}), 0.0, node, true);
  check(Gc({0, 2}, {0, 3}), gc_ref_2);

  boost::multi::array<Type, 1, shared_allocator<Type>> SMbuff2(iextensions<1u>{NMO * (NMO + NEL)},
                                                               shared_allocator<Type>{node_});

  array_ref G2(to_address(SMbuff2.origin()), {NMO, NMO});
  array_ref Gc2(to_address(SMbuff2.origin()) + NMO * NMO, {NEL, NMO});

  // switch comm

  //MAM
  // turning off this test, somehow failing in rhea with gcc5.3
  // not sure why!!! Will fix soon!
  // This is probably a bug in the native implementation of csrmm!!!
  /*
  ov_=SDet.MixedDensityMatrix(Acsr,B,G2,0.0,node_,false); 
  check(G2,g_ref);
  ov_=SDet.MixedDensityMatrix(Acsr[{0,2,0,3}],
                          B({0,3},{0,2}),
                          G2({0,3},{0,3}),0.0,node_,false);
  check(G2({0,3},{0,3}),g_ref_2);
  ov_=SDet.MixedDensityMatrix(Acsr,B,Gc2,0.0,node_,true); check(Gc2,gc_ref);
  ov_=SDet.MixedDensityMatrix(Acsr[{0,2,0,3}],
                          B({0,3},{0,2}),
                          Gc2({0,2},{0,3}),0.0,node_,true);
  check(Gc2({0,2},{0,3}),gc_ref_2);
*/
  release_memory_managers();
}


TEST_CASE("SDetOps_complex_serial", "[sdet_ops]")
{
  auto world = boost::mpi3::environment::get_world_instance();
  auto node  = world.split_shared(world.rank());

#if defined(ENABLE_CUDA) || defined(ENABLE_HIP)
  arch::INIT(node);
  using Alloc = device::device_allocator<ComplexType>;
#else
  using Alloc = std::allocator<ComplexType>;
#endif
  setup_memory_managers(node, 10uL * 1024uL * 1024uL);
  SDetOps_complex_serial<Alloc, DeviceBufferManager>(Alloc{}, DeviceBufferManager{});
  release_memory_managers();
}

/*
TEST_CASE("SDetOps_complex_mpi3", "[sdet_ops]")
{


}

TEST_CASE("SDetOps_complex_csr", "[sdet_ops]")
{

}
*/
} // namespace afqmc
} // namespace qmcplusplus
