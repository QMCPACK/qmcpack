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


#include "Message/catch_mpi_main.hpp"

//#define CATCH_CONFIG_MAIN
//#include "catch.hpp"
#include "Configuration.h"

#include "OhmmsApp/ProjectData.h"

// Avoid the need to link with other libraries just to get APP_ABORT
#undef APP_ABORT
#define APP_ABORT(x) std::cout << x; exit(0);

#include <stdio.h>
#include <string>
#include <vector>
#include <complex>

#undef NDEBUG

#include "mpi3/shared_communicator.hpp"
#include "mpi3/environment.hpp"

#include "boost/multi_array.hpp"
#include "AFQMC/SlaterDeterminantOperations/SlaterDetOperations.hpp"
#include "AFQMC/SlaterDeterminantOperations/mixed_density_matrix.hpp"

using std::string;
using std::complex;
using std::cout;
using std::endl;

using boost::extents;
using boost::indices;
using range_t = boost::multi_array_types::index_range;

namespace qmcplusplus
{

template<class M1, class M2>
void check(M1&& A, M2& B)
{
  REQUIRE(A.shape()[0] == B.shape()[0]);
  REQUIRE(A.shape()[1] == B.shape()[1]);
  for(int i=0; i<A.shape()[0]; i++)
    for(int j=0; j<A.shape()[1]; j++)
      REQUIRE(A[i][j] == Approx(B[i][j]));
}

using namespace afqmc;

TEST_CASE("SDetOps_double_serial", "[sdet_ops]")
{
  Communicate *c;
  OHMMS::Controller->initialize(0, NULL);
  //c = OHMMS::Controller;

  const int NMO = 4;
  const int NEL = 3;

  using Type = double;
  using vector = std::vector<Type>;
  using multi_array = boost::multi_array<Type,2>;
  using multi_array_ref = boost::multi_array_ref<Type,2>;

  // some arbitrary matrices
  vector m_a = {
   0.90000,   2.40000,   3.00000,
   0.40000,   1.00000,   1.20000,
   1.40000,   1.60000,   3.60000,
   0.40000,   0.20000,   0.10000
  };
  vector m_b = {
   1.90000,   1.40000,   0.40000,
   1.40000,   0.20000,   2.20000,
   0.40000,   2.60000,   0.60000,
   1.10000,   0.30000,   0.90000
  };

  multi_array A(extents[NMO][NEL]);
  multi_array B(extents[NMO][NEL]);

  for(int i=0, k=0; i<A.shape()[0]; i++)
    for(int j=0; j<A.shape()[1]; j++,k++)
       A[i][j] = m_a[k];

  for(int i=0, k=0; i<B.shape()[0]; i++)
    for(int j=0; j<B.shape()[1]; j++,k++)
       B[i][j] = m_b[k];

  multi_array_ref Aref(m_a.data(),extents[NMO][NEL]);
  multi_array_ref Bref(m_b.data(),extents[NMO][NEL]);

  SlaterDetOperations<Type> SDet(NMO,NEL);

  /**** Overlaps ****/
  REQUIRE(SDet.Overlap(A,B) == Approx(5.10443199999999));
  REQUIRE(SDet.Overlap(Aref,B) == Approx(5.10443199999999));
  REQUIRE(SDet.Overlap(A,Bref) == Approx(5.10443199999999));
  REQUIRE(SDet.Overlap(Aref,Bref) == Approx(5.10443199999999));

  // Test array_view
  REQUIRE(SDet.Overlap(A[indices[range_t()][range_t()]],B) == Approx(5.10443199999999));
  REQUIRE(SDet.Overlap(A,B[indices[range_t()][range_t()]]) == Approx(5.10443199999999));

  multi_array A_ = A[indices[range_t(0,3)][range_t(0,2)]];
  multi_array B_ = B[indices[range_t(0,3)][range_t(0,2)]];
  REQUIRE(SDet.Overlap(A[indices[range_t(0,3)][range_t(0,2)]],
                       B[indices[range_t(0,3)][range_t(0,2)]]) == Approx(-11.0204000000000));
  REQUIRE(SDet.Overlap(A[indices[range_t(0,3)][range_t(0,2)]],B_) == Approx(-11.0204000000000));
  REQUIRE(SDet.Overlap(A_,B[indices[range_t(0,3)][range_t(0,2)]]) == Approx(-11.0204000000000));


  /**** Density Matrices *****/
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

  multi_array_ref g_ref(v_ref.data(),extents[NMO][NMO]);
  multi_array_ref gc_ref(vc_ref.data(),extents[NEL][NMO]);
  multi_array_ref g_ref_2(v_ref_2.data(),extents[3][3]);
  multi_array_ref gc_ref_2(vc_ref_2.data(),extents[2][3]);

  multi_array G(extents[NMO][NMO]);
  multi_array Gc(extents[NEL][NMO]);

  SDet.MixedDensityMatrix(A,B,G,false); check(G,g_ref);
  SDet.MixedDensityMatrix(Aref,B,G,false); check(G,g_ref);
  SDet.MixedDensityMatrix(A,Bref,G,false); check(G,g_ref);
  SDet.MixedDensityMatrix(Aref,Bref,G,false); check(G,g_ref);

  SDet.MixedDensityMatrix(A[indices[range_t(0,3)][range_t(0,2)]],
                          B[indices[range_t(0,3)][range_t(0,2)]],
                          G[indices[range_t(0,3)][range_t(0,3)]],false);
  check(G[indices[range_t(0,3)][range_t(0,3)]],g_ref_2);
  SDet.MixedDensityMatrix(A_,
                          B[indices[range_t(0,3)][range_t(0,2)]],
                          G[indices[range_t(0,3)][range_t(0,3)]],false);
  check(G[indices[range_t(0,3)][range_t(0,3)]],g_ref_2);
  SDet.MixedDensityMatrix(A[indices[range_t(0,3)][range_t(0,2)]],
                          B_,
                          G[indices[range_t(0,3)][range_t(0,3)]],false);
  check(G[indices[range_t(0,3)][range_t(0,3)]],g_ref_2);

  SDet.MixedDensityMatrix(A,B,Gc,true); check(Gc,gc_ref);
  SDet.MixedDensityMatrix(Aref,B,Gc,true); check(Gc,gc_ref);
  SDet.MixedDensityMatrix(A,Bref,Gc,true); check(Gc,gc_ref);
  SDet.MixedDensityMatrix(Aref,Bref,Gc,true); check(Gc,gc_ref);

  SDet.MixedDensityMatrix(A[indices[range_t(0,3)][range_t(0,2)]],
                          B[indices[range_t(0,3)][range_t(0,2)]],
                          Gc[indices[range_t(0,2)][range_t(0,3)]],true);
  check(Gc[indices[range_t(0,2)][range_t(0,3)]],gc_ref_2);
  SDet.MixedDensityMatrix(A_,
                          B[indices[range_t(0,3)][range_t(0,2)]],
                          Gc[indices[range_t(0,2)][range_t(0,3)]],true);
  check(Gc[indices[range_t(0,2)][range_t(0,3)]],gc_ref_2);
  SDet.MixedDensityMatrix(A[indices[range_t(0,3)][range_t(0,2)]],
                          B_,
                          Gc[indices[range_t(0,2)][range_t(0,3)]],true);
  check(Gc[indices[range_t(0,2)][range_t(0,3)]],gc_ref_2);

  multi_array Q=A;

  // Orthogonalize
  Type detR = SDet.Orthogonalize(Q);
  REQUIRE( SDet.Overlap(Q,Q) == Approx(1.0)  );

}

TEST_CASE("SDetOps_double_mpi3", "[sdet_ops]")
{

  Communicate *c = OHMMS::Controller;

  using boost::mpi3::shared_communicator;

  shared_communicator node = c->comm.split_shared(c->comm.rank());

  const int NMO = 4;
  const int NEL = 3;

  using Type = double;
  using vector = std::vector<Type>;
  using multi_array = boost::multi_array<Type,2>;
  using multi_array_ref = boost::multi_array_ref<Type,2>;

  const Type ov = 5.10443199999999;
  const Type ov2 = -11.0204000000000;


  // some arbitrary matrices
  vector m_a = {
   0.90000,   2.40000,   3.00000,
   0.40000,   1.00000,   1.20000,
   1.40000,   1.60000,   3.60000,
   0.40000,   0.20000,   0.10000
  };
  vector m_b = {
   1.90000,   1.40000,   0.40000,
   1.40000,   0.20000,   2.20000,
   0.40000,   2.60000,   0.60000,
   1.10000,   0.30000,   0.90000
  };

  multi_array A(extents[NMO][NEL]);
  multi_array B(extents[NMO][NEL]);

  for(int i=0, k=0; i<A.shape()[0]; i++)
    for(int j=0; j<A.shape()[1]; j++,k++)
       A[i][j] = m_a[k];

  for(int i=0, k=0; i<B.shape()[0]; i++)
    for(int j=0; j<B.shape()[1]; j++,k++)
       B[i][j] = m_b[k];

  multi_array_ref Aref(m_a.data(),extents[NMO][NEL]);
  multi_array_ref Bref(m_b.data(),extents[NMO][NEL]);

  SlaterDetOperations<Type> SDet(NMO,NEL);

  /**** Overlaps ****/
  REQUIRE(SDet.Overlap(A,B,node) == Approx(5.10443199999999));
  REQUIRE(SDet.Overlap(Aref,B,node) == Approx(5.10443199999999));
  REQUIRE(SDet.Overlap(A,Bref,node) == Approx(5.10443199999999));
  REQUIRE(SDet.Overlap(Aref,Bref,node) == Approx(5.10443199999999));

  // Test array_view
  REQUIRE(SDet.Overlap(A[indices[range_t()][range_t()]],B,node) == Approx(5.10443199999999));
  REQUIRE(SDet.Overlap(A,B[indices[range_t()][range_t()]],node) == Approx(5.10443199999999));

  multi_array A_ = A[indices[range_t(0,3)][range_t(0,2)]];
  multi_array B_ = B[indices[range_t(0,3)][range_t(0,2)]];
  REQUIRE(SDet.Overlap(A[indices[range_t(0,3)][range_t(0,2)]],
                       B[indices[range_t(0,3)][range_t(0,2)]],node) == Approx(-11.0204000000000));
  REQUIRE(SDet.Overlap(A[indices[range_t(0,3)][range_t(0,2)]],B_) == Approx(-11.0204000000000));
  REQUIRE(SDet.Overlap(A_,B[indices[range_t(0,3)][range_t(0,2)]],node) == Approx(-11.0204000000000));

  shared_communicator node_ = node.split(node.rank()%2);
  REQUIRE(SDet.Overlap(A,B,node_) == Approx(5.10443199999999));
  REQUIRE(SDet.Overlap(A[indices[range_t(0,3)][range_t(0,2)]],
                       B[indices[range_t(0,3)][range_t(0,2)]],node_) == Approx(-11.0204000000000));


  /**** Density Matrices *****/
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

  multi_array_ref g_ref(v_ref.data(),extents[NMO][NMO]);
  multi_array_ref gc_ref(vc_ref.data(),extents[NEL][NMO]);
  multi_array_ref g_ref_2(v_ref_2.data(),extents[3][3]);
  multi_array_ref gc_ref_2(vc_ref_2.data(),extents[2][3]);

  using SHM_Buffer = mpi3_SHMBuffer<Type>;
  SHM_Buffer SMbuff(node,NMO*(NMO+NEL));

  multi_array_ref G(SMbuff.data(),extents[NMO][NMO]);
  multi_array_ref Gc(SMbuff.data()+NMO*NMO,extents[NEL][NMO]);

  SDet.MixedDensityMatrix(A,B,G,node,false); check(G,g_ref);
  SDet.MixedDensityMatrix(Aref,B,G,node,false); check(G,g_ref);
  SDet.MixedDensityMatrix(A,Bref,G,node,false); check(G,g_ref);
  SDet.MixedDensityMatrix(Aref,Bref,G,node,false); check(G,g_ref);

  SDet.MixedDensityMatrix(A[indices[range_t(0,3)][range_t(0,2)]],
                          B[indices[range_t(0,3)][range_t(0,2)]],
                          G[indices[range_t(0,3)][range_t(0,3)]],node,false);
  check(G[indices[range_t(0,3)][range_t(0,3)]],g_ref_2);
  SDet.MixedDensityMatrix(A_,
                          B[indices[range_t(0,3)][range_t(0,2)]],
                          G[indices[range_t(0,3)][range_t(0,3)]],node,false);
  check(G[indices[range_t(0,3)][range_t(0,3)]],g_ref_2);
  SDet.MixedDensityMatrix(A[indices[range_t(0,3)][range_t(0,2)]],
                          B_,
                          G[indices[range_t(0,3)][range_t(0,3)]],node,false);
  check(G[indices[range_t(0,3)][range_t(0,3)]],g_ref_2);

  SDet.MixedDensityMatrix(A,B,Gc,node,true); check(Gc,gc_ref);
  SDet.MixedDensityMatrix(Aref,B,Gc,node,true); check(Gc,gc_ref);
  SDet.MixedDensityMatrix(A,Bref,Gc,node,true); check(Gc,gc_ref);
  SDet.MixedDensityMatrix(Aref,Bref,Gc,node,true); check(Gc,gc_ref);

  SDet.MixedDensityMatrix(A[indices[range_t(0,3)][range_t(0,2)]],
                          B[indices[range_t(0,3)][range_t(0,2)]],
                          Gc[indices[range_t(0,2)][range_t(0,3)]],node,true);
  check(Gc[indices[range_t(0,2)][range_t(0,3)]],gc_ref_2);
  SDet.MixedDensityMatrix(A_,
                          B[indices[range_t(0,3)][range_t(0,2)]],
                          Gc[indices[range_t(0,2)][range_t(0,3)]],node,true);
  check(Gc[indices[range_t(0,2)][range_t(0,3)]],gc_ref_2);
  SDet.MixedDensityMatrix(A[indices[range_t(0,3)][range_t(0,2)]],
                          B_,
                          Gc[indices[range_t(0,2)][range_t(0,3)]],node,true);
  check(Gc[indices[range_t(0,2)][range_t(0,3)]],gc_ref_2);

  SHM_Buffer SMbuff2(node_,NMO*(NMO+NEL));

  multi_array_ref G2(SMbuff2.data(),extents[NMO][NMO]);
  multi_array_ref Gc2(SMbuff2.data()+NMO*NMO,extents[NEL][NMO]);

  // switch comm
  SDet.MixedDensityMatrix(A,B,G2,node_,false); check(G2,g_ref);
  SDet.MixedDensityMatrix(A[indices[range_t(0,3)][range_t(0,2)]],
                          B[indices[range_t(0,3)][range_t(0,2)]],
                          G2[indices[range_t(0,3)][range_t(0,3)]],node_,false);
  check(G2[indices[range_t(0,3)][range_t(0,3)]],g_ref_2);
  SDet.MixedDensityMatrix(A,B,Gc2,node_,true); check(Gc2,gc_ref);
  SDet.MixedDensityMatrix(A[indices[range_t(0,3)][range_t(0,2)]],
                          B[indices[range_t(0,3)][range_t(0,2)]],
                          Gc2[indices[range_t(0,2)][range_t(0,3)]],node_,true);
  check(Gc2[indices[range_t(0,2)][range_t(0,3)]],gc_ref_2);

}

}
