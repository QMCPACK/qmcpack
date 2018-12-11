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

#include "Message/catch_mpi_main.hpp"

//#define CATCH_CONFIG_MAIN
//#include "catch.hpp"
#include "Configuration.h"

#include "OhmmsApp/ProjectData.h"

// Avoid the need to link with other libraries just to get APP_ABORT
#undef APP_ABORT
#define APP_ABORT(x) {std::cout << x <<std::endl; exit(0);}

#include <stdio.h>
#include <string>
#include <vector>
#include <complex>

#include "mpi3/shared_communicator.hpp"
#include "mpi3/environment.hpp"

#include "boost/multi_array.hpp"
#include "AFQMC/Matrix/csr_matrix_construct.hpp"
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

void myREQUIRE(const double& a, const double& b)
{
  REQUIRE(a == Approx(b));
}

void myREQUIRE(const std::complex<double>& a, const std::complex<double>& b)
{
  REQUIRE(a.real() == Approx(b.real()));
  REQUIRE(a.imag() == Approx(b.imag()));
}

template<class M1, class M2>
void check(M1&& A, M2& B)
{
  REQUIRE(A.shape()[0] == B.shape()[0]);
  REQUIRE(A.shape()[1] == B.shape()[1]);
  for(int i=0; i<A.shape()[0]; i++)
    for(int j=0; j<A.shape()[1]; j++)
      myREQUIRE(A[i][j],B[i][j]);
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

  multi_array A(extents[NEL][NMO]);
  multi_array B(extents[NMO][NEL]);

  for(int i=0, k=0; i<A.shape()[0]; i++)
    for(int j=0; j<A.shape()[1]; j++,k++)
       A[i][j] = m_a[k]; 

  for(int i=0, k=0; i<B.shape()[0]; i++)
    for(int j=0; j<B.shape()[1]; j++,k++)
       B[i][j] = m_b[k];

  multi_array_ref Aref(m_a.data(),extents[NEL][NMO]);
  multi_array_ref Bref(m_b.data(),extents[NMO][NEL]);

  SlaterDetOperations<Type> SDet(NMO,NEL);

  /**** Overlaps ****/  
  REQUIRE(SDet.Overlap(A,B) == Approx(ov)); 
  REQUIRE(SDet.Overlap(Aref,B) == Approx(ov)); 
  REQUIRE(SDet.Overlap(A,Bref) == Approx(ov)); 
  REQUIRE(SDet.Overlap(Aref,Bref) == Approx(ov)); 

  // Test array_view
  REQUIRE(SDet.Overlap(A[indices[range_t()][range_t()]],B) == Approx(ov));
  REQUIRE(SDet.Overlap(A,B[indices[range_t()][range_t()]]) == Approx(ov));

  multi_array A_ = A[indices[range_t(0,2)][range_t(0,3)]];
  multi_array B_ = B[indices[range_t(0,3)][range_t(0,2)]];
  REQUIRE(SDet.Overlap(A[indices[range_t(0,2)][range_t(0,3)]],
                       B[indices[range_t(0,3)][range_t(0,2)]]) == Approx(ov2));
  REQUIRE(SDet.Overlap(A[indices[range_t(0,2)][range_t(0,3)]],B_) == Approx(ov2));
  REQUIRE(SDet.Overlap(A_,B[indices[range_t(0,3)][range_t(0,2)]]) == Approx(ov2));


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

  SDet.MixedDensityMatrix(A[indices[range_t(0,2)][range_t(0,3)]],
                          B[indices[range_t(0,3)][range_t(0,2)]],
                          G[indices[range_t(0,3)][range_t(0,3)]],false); 
  check(G[indices[range_t(0,3)][range_t(0,3)]],g_ref_2);  
  SDet.MixedDensityMatrix(A_,
                          B[indices[range_t(0,3)][range_t(0,2)]],
                          G[indices[range_t(0,3)][range_t(0,3)]],false);
  check(G[indices[range_t(0,3)][range_t(0,3)]],g_ref_2);
  SDet.MixedDensityMatrix(A[indices[range_t(0,2)][range_t(0,3)]],
                          B_,
                          G[indices[range_t(0,3)][range_t(0,3)]],false);
  check(G[indices[range_t(0,3)][range_t(0,3)]],g_ref_2);

  SDet.MixedDensityMatrix(A,B,Gc,true); check(Gc,gc_ref);
  SDet.MixedDensityMatrix(Aref,B,Gc,true); check(Gc,gc_ref);
  SDet.MixedDensityMatrix(A,Bref,Gc,true); check(Gc,gc_ref);
  SDet.MixedDensityMatrix(Aref,Bref,Gc,true); check(Gc,gc_ref);  

  SDet.MixedDensityMatrix(A[indices[range_t(0,2)][range_t(0,3)]],
                          B[indices[range_t(0,3)][range_t(0,2)]],
                          Gc[indices[range_t(0,2)][range_t(0,3)]],true);
  check(Gc[indices[range_t(0,2)][range_t(0,3)]],gc_ref_2);
  SDet.MixedDensityMatrix(A_,
                          B[indices[range_t(0,3)][range_t(0,2)]],
                          Gc[indices[range_t(0,2)][range_t(0,3)]],true);
  check(Gc[indices[range_t(0,2)][range_t(0,3)]],gc_ref_2);
  SDet.MixedDensityMatrix(A[indices[range_t(0,2)][range_t(0,3)]],
                          B_,
                          Gc[indices[range_t(0,2)][range_t(0,3)]],true);
  check(Gc[indices[range_t(0,2)][range_t(0,3)]],gc_ref_2);

  multi_array Q=B;

  // Orthogonalize
  Type detR = SDet.Orthogonalize(Q);
  REQUIRE( SDet.Overlap_noHerm(Q,Q) == Approx(1.0)  );

}

TEST_CASE("SDetOps_double_mpi3", "[sdet_ops]")
{

  Communicate *c;
  OHMMS::Controller->initialize(0, NULL);
  //c = OHMMS::Controller;
  
  using boost::mpi3::shared_communicator;

  shared_communicator node = boost::mpi3::world.split_shared(boost::mpi3::world.rank()); 

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

  multi_array A(extents[NEL][NMO]);
  multi_array B(extents[NMO][NEL]);

  for(int i=0, k=0; i<A.shape()[0]; i++)
    for(int j=0; j<A.shape()[1]; j++,k++)
       A[i][j] = m_a[k];

  for(int i=0, k=0; i<B.shape()[0]; i++)
    for(int j=0; j<B.shape()[1]; j++,k++)
       B[i][j] = m_b[k];

  multi_array_ref Aref(m_a.data(),extents[NEL][NMO]);
  multi_array_ref Bref(m_b.data(),extents[NMO][NEL]);

  SlaterDetOperations<Type> SDet(NMO,NEL);

  /**** Overlaps ****/
  REQUIRE(SDet.Overlap(A,B,node) == Approx(ov));
  REQUIRE(SDet.Overlap(Aref,B,node) == Approx(ov));
  REQUIRE(SDet.Overlap(A,Bref,node) == Approx(ov));
  REQUIRE(SDet.Overlap(Aref,Bref,node) == Approx(ov));

  // Test array_view
  REQUIRE(SDet.Overlap(A[indices[range_t()][range_t()]],B,node) == Approx(ov));
  REQUIRE(SDet.Overlap(A,B[indices[range_t()][range_t()]],node) == Approx(ov));

  multi_array A_ = A[indices[range_t(0,2)][range_t(0,3)]];
  multi_array B_ = B[indices[range_t(0,3)][range_t(0,2)]];
  REQUIRE(SDet.Overlap(A[indices[range_t(0,2)][range_t(0,3)]],
                       B[indices[range_t(0,3)][range_t(0,2)]],node) == Approx(ov2));
  REQUIRE(SDet.Overlap(A[indices[range_t(0,2)][range_t(0,3)]],B_) == Approx(ov2));
  REQUIRE(SDet.Overlap(A_,B[indices[range_t(0,3)][range_t(0,2)]],node) == Approx(ov2));

  shared_communicator node_ = node.split(node.rank()%2); 
  REQUIRE(SDet.Overlap(A,B,node_) == Approx(ov));
  REQUIRE(SDet.Overlap(A[indices[range_t(0,2)][range_t(0,3)]],
                       B[indices[range_t(0,3)][range_t(0,2)]],node_) == Approx(ov2));


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

  SDet.MixedDensityMatrix(A[indices[range_t(0,2)][range_t(0,3)]],
                          B[indices[range_t(0,3)][range_t(0,2)]],
                          G[indices[range_t(0,3)][range_t(0,3)]],node,false);
  check(G[indices[range_t(0,3)][range_t(0,3)]],g_ref_2);
  SDet.MixedDensityMatrix(A_,
                          B[indices[range_t(0,3)][range_t(0,2)]],
                          G[indices[range_t(0,3)][range_t(0,3)]],node,false);
  check(G[indices[range_t(0,3)][range_t(0,3)]],g_ref_2);
  SDet.MixedDensityMatrix(A[indices[range_t(0,2)][range_t(0,3)]],
                          B_,
                          G[indices[range_t(0,3)][range_t(0,3)]],node,false);
  check(G[indices[range_t(0,3)][range_t(0,3)]],g_ref_2);

  SDet.MixedDensityMatrix(A,B,Gc,node,true); check(Gc,gc_ref);
  SDet.MixedDensityMatrix(Aref,B,Gc,node,true); check(Gc,gc_ref);
  SDet.MixedDensityMatrix(A,Bref,Gc,node,true); check(Gc,gc_ref);
  SDet.MixedDensityMatrix(Aref,Bref,Gc,node,true); check(Gc,gc_ref);

  SDet.MixedDensityMatrix(A[indices[range_t(0,2)][range_t(0,3)]],
                          B[indices[range_t(0,3)][range_t(0,2)]],
                          Gc[indices[range_t(0,2)][range_t(0,3)]],node,true);
  check(Gc[indices[range_t(0,2)][range_t(0,3)]],gc_ref_2);
  SDet.MixedDensityMatrix(A_,
                          B[indices[range_t(0,3)][range_t(0,2)]],
                          Gc[indices[range_t(0,2)][range_t(0,3)]],node,true);
  check(Gc[indices[range_t(0,2)][range_t(0,3)]],gc_ref_2);
  SDet.MixedDensityMatrix(A[indices[range_t(0,2)][range_t(0,3)]],
                          B_,
                          Gc[indices[range_t(0,2)][range_t(0,3)]],node,true);
  check(Gc[indices[range_t(0,2)][range_t(0,3)]],gc_ref_2);

  SHM_Buffer SMbuff2(node_,NMO*(NMO+NEL));

  multi_array_ref G2(SMbuff2.data(),extents[NMO][NMO]);
  multi_array_ref Gc2(SMbuff2.data()+NMO*NMO,extents[NEL][NMO]);

  // switch comm
  SDet.MixedDensityMatrix(A,B,G2,node_,false); check(G2,g_ref);
  SDet.MixedDensityMatrix(A[indices[range_t(0,2)][range_t(0,3)]],
                          B[indices[range_t(0,3)][range_t(0,2)]],
                          G2[indices[range_t(0,3)][range_t(0,3)]],node_,false);
  check(G2[indices[range_t(0,3)][range_t(0,3)]],g_ref_2);
  SDet.MixedDensityMatrix(A,B,Gc2,node_,true); check(Gc2,gc_ref);
  SDet.MixedDensityMatrix(A[indices[range_t(0,2)][range_t(0,3)]],
                          B[indices[range_t(0,3)][range_t(0,2)]],
                          Gc2[indices[range_t(0,2)][range_t(0,3)]],node_,true);
  check(Gc2[indices[range_t(0,2)][range_t(0,3)]],gc_ref_2);  
   
}

TEST_CASE("SDetOps_complex_serial", "[sdet_ops]")
{
  Communicate *c;
  OHMMS::Controller->initialize(0, NULL);
  //c = OHMMS::Controller;

  const int NMO = 4;
  const int NEL = 3;  

  using Type = std::complex<double>;
  using vector = std::vector<Type>;
  using multi_array = boost::multi_array<Type,2>;
  using multi_array_ref = boost::multi_array_ref<Type,2>;
  using namespace std::complex_literals;

  const Type ov = -7.62332599999999 + 22.20453200000000i; 
  const Type ov2 = -10.37150000000000 -  7.15750000000000i; 

  // some arbitrary matrices
  vector m_a = {
//   0.90000 + 0.10000i,   2.40000 + 0.20000i,   3.00000 + 0.30000i,
//   0.40000 + 0.40000i,   1.00000 + 0.50000i,   1.20000 + 0.10000i,
//   1.40000 + 0.20000i,   1.60000 + 0.30000i,   3.60000 + 0.40000i,
//   0.40000 + 0.50000i,   0.20000 + 0.10000i,   0.10000 + 0.20000i
    0.90000 +0.10000i , 0.40000 + 0.40000i, 1.40000 + 0.20000i , 0.40000 + 0.50000i,
    2.40000 +0.20000i , 1.00000 + 0.50000i, 1.60000 + 0.30000i , 0.20000 + 0.10000i,
    3.00000 +0.30000i , 1.20000 + 0.10000i, 3.60000 + 0.40000i , 0.10000 + 0.20000i
  };  
  vector m_b = {
   1.90000 + 0.60000i,   1.40000 + 0.70000i,   0.40000 + 0.80000i,
   1.40000 + 0.90000i,   0.20000 + 0.50000i,   2.20000 + 0.60000i,
   0.40000 + 0.70000i,   2.60000 + 0.80000i,   0.60000 + 0.90000i,
   1.10000 + 0.50000i,   0.30000 + 0.60000i,   0.90000 + 0.70000i
  };

  multi_array A(extents[NEL][NMO]);
  multi_array B(extents[NMO][NEL]);

  for(int i=0, k=0; i<A.shape()[0]; i++)
    for(int j=0; j<A.shape()[1]; j++,k++)
       A[i][j] = m_a[k]; 

  for(int i=0, k=0; i<B.shape()[0]; i++)
    for(int j=0; j<B.shape()[1]; j++,k++)
       B[i][j] = m_b[k];

  multi_array_ref Aref(m_a.data(),extents[NEL][NMO]);
  multi_array_ref Bref(m_b.data(),extents[NMO][NEL]);

  SlaterDetOperations<Type> SDet(NMO,NEL);

  /**** Overlaps ****/  
  myREQUIRE(SDet.Overlap(A,B),ov); 
  myREQUIRE(SDet.Overlap(Aref,B),ov);
  myREQUIRE(SDet.Overlap(A,Bref),ov);
  myREQUIRE(SDet.Overlap(Aref,Bref),ov);

  // Test array_view
  myREQUIRE(SDet.Overlap(A[indices[range_t()][range_t()]],B),ov);
  myREQUIRE(SDet.Overlap(A,B[indices[range_t()][range_t()]]),ov);

  multi_array A_ = A[indices[range_t(0,2)][range_t(0,3)]];
  multi_array B_ = B[indices[range_t(0,3)][range_t(0,2)]];
  myREQUIRE(SDet.Overlap(A[indices[range_t(0,2)][range_t(0,3)]],
                       B[indices[range_t(0,3)][range_t(0,2)]]),ov2);
  myREQUIRE(SDet.Overlap(A[indices[range_t(0,2)][range_t(0,3)]],B_),ov2);
  myREQUIRE(SDet.Overlap(A_,B[indices[range_t(0,3)][range_t(0,2)]]),ov2);

  /**** Density Matrices *****/
  vector v_ref = {
   1.17573619385025996 - 0.01580426445014660i,  -0.25295981756593167 + 0.28594469607401085i,
  -0.07724502823533341 - 0.09687959052155870i,   0.30512858581808422 - 0.04506898328729603i,

   0.17912592806889663 + 0.08374315906672802i,   0.59381451118048767 + 0.13438888951771200i,
  -0.02021475320201610 - 0.13737561982083193i,   0.32095003919745313 + 0.12832750636154097i,

  -0.04919549564646425 + 0.05402065222741825i,  -0.00286990878355775 - 0.15806420733175885i,
   1.05069081188635494 + 0.00793429988912356i,  -0.08048239150997794 + 0.09917405634760490i,

  -0.46434598219794548 - 0.09896890422706731i,   0.87746748807670427 - 0.53417787485950319i,
   0.12162735438647077 + 0.31042401735800573i,   0.17975848308289613 - 0.12651892495668898i
  };
  vector vc_ref = {
  -0.8721879971495297 + 0.6787593377585239i,   0.3244278932250768 - 1.8083537898275881i,
   0.6130713546530860 + 0.0399736955598931i,   0.0132562806444336 - 0.2882495766584950i,

   0.8853626557603183 + 0.0978868569224204i,  -0.0598704127345155 - 0.0470889603064014i,
  -0.7392693424168300 - 0.0715317395994149i,   0.3721269544963505 - 0.1797896522886788i,

  -0.0567190307022984 - 0.3114847157576828i,  -0.1290126440468128 + 0.6815705660308808i,
   0.3787012855335005 + 0.0039188686237135i,  -0.2005543456941538 + 0.2100886953142371i
  };  
  vector v_ref_2 = {
   0.7361983496013835 - 0.0956505507662245i,   0.6467449689807925 + 0.2297471806893873i,
   0.0189270005620390 - 0.1727975708935829i,

   0.2986893588843604 - 0.0099955730815030i,   0.2696068479051557 + 0.0292039710860386i,
   0.0497734835066391 + 0.1783200397050796i,

   0.1020030246826092 + 0.0344707468383766i,  -0.2500340021988402 - 0.0826863644855427i,
   0.9941948024934623 + 0.0664465796801866i
  };  
  vector vc_ref_2 = {
  -0.489369975192701 + 0.103038673040713i,  -0.858850485405126 - 0.275734238124941i,
   1.219791948842170 + 0.118626922447301i,
   0.486337033653898 - 0.098631569047656i,   0.595497821653010 + 0.185288949671560i,
  -0.455373025165330 - 0.129360996228044i
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

  SDet.MixedDensityMatrix(A[indices[range_t(0,2)][range_t(0,3)]],
                          B[indices[range_t(0,3)][range_t(0,2)]],
                          G[indices[range_t(0,3)][range_t(0,3)]],false); 
  check(G[indices[range_t(0,3)][range_t(0,3)]],g_ref_2);  
  SDet.MixedDensityMatrix(A_,
                          B[indices[range_t(0,3)][range_t(0,2)]],
                          G[indices[range_t(0,3)][range_t(0,3)]],false);
  check(G[indices[range_t(0,3)][range_t(0,3)]],g_ref_2);
  SDet.MixedDensityMatrix(A[indices[range_t(0,2)][range_t(0,3)]],
                          B_,
                          G[indices[range_t(0,3)][range_t(0,3)]],false);
  check(G[indices[range_t(0,3)][range_t(0,3)]],g_ref_2);

  SDet.MixedDensityMatrix(A,B,Gc,true); check(Gc,gc_ref);
  SDet.MixedDensityMatrix(Aref,B,Gc,true); check(Gc,gc_ref);
  SDet.MixedDensityMatrix(A,Bref,Gc,true); check(Gc,gc_ref);
  SDet.MixedDensityMatrix(Aref,Bref,Gc,true); check(Gc,gc_ref);  

  SDet.MixedDensityMatrix(A[indices[range_t(0,2)][range_t(0,3)]],
                          B[indices[range_t(0,3)][range_t(0,2)]],
                          Gc[indices[range_t(0,2)][range_t(0,3)]],true);
  check(Gc[indices[range_t(0,2)][range_t(0,3)]],gc_ref_2);
  SDet.MixedDensityMatrix(A_,
                          B[indices[range_t(0,3)][range_t(0,2)]],
                          Gc[indices[range_t(0,2)][range_t(0,3)]],true);
  check(Gc[indices[range_t(0,2)][range_t(0,3)]],gc_ref_2);
  SDet.MixedDensityMatrix(A[indices[range_t(0,2)][range_t(0,3)]],
                          B_,
                          Gc[indices[range_t(0,2)][range_t(0,3)]],true);
  check(Gc[indices[range_t(0,2)][range_t(0,3)]],gc_ref_2);

  multi_array Q=B;

  // Orthogonalize
  Type detR = SDet.Orthogonalize(Q);
  myREQUIRE( SDet.Overlap_noHerm(Q,Q), std::complex<double>(1.,0.));

}

TEST_CASE("SDetOps_complex_mpi3", "[sdet_ops]")
{

  Communicate *c;
  OHMMS::Controller->initialize(0, NULL);
  //c = OHMMS::Controller;
  
  using boost::mpi3::shared_communicator;

  shared_communicator node = boost::mpi3::world.split_shared(boost::mpi3::world.rank()); 

  const int NMO = 4;
  const int NEL = 3;

  using Type = std::complex<double>;
  using vector = std::vector<Type>;
  using multi_array = boost::multi_array<Type,2>;
  using multi_array_ref = boost::multi_array_ref<Type,2>;
  using namespace std::complex_literals;

  const Type ov = -7.62332599999999 + 22.20453200000000i;
  const Type ov2 = -10.37150000000000 -  7.15750000000000i;

  // some arbitrary matrices
  vector m_a = {
//   0.90000 + 0.10000i,   2.40000 + 0.20000i,   3.00000 + 0.30000i,
//   0.40000 + 0.40000i,   1.00000 + 0.50000i,   1.20000 + 0.10000i,
//   1.40000 + 0.20000i,   1.60000 + 0.30000i,   3.60000 + 0.40000i,
//   0.40000 + 0.50000i,   0.20000 + 0.10000i,   0.10000 + 0.20000i
    0.90000 +0.10000i , 0.40000 + 0.40000i, 1.40000 + 0.20000i , 0.40000 + 0.50000i,
    2.40000 +0.20000i , 1.00000 + 0.50000i, 1.60000 + 0.30000i , 0.20000 + 0.10000i,
    3.00000 +0.30000i , 1.20000 + 0.10000i, 3.60000 + 0.40000i , 0.10000 + 0.20000i
  };
  vector m_b = {
   1.90000 + 0.60000i,   1.40000 + 0.70000i,   0.40000 + 0.80000i,
   1.40000 + 0.90000i,   0.20000 + 0.50000i,   2.20000 + 0.60000i,
   0.40000 + 0.70000i,   2.60000 + 0.80000i,   0.60000 + 0.90000i,
   1.10000 + 0.50000i,   0.30000 + 0.60000i,   0.90000 + 0.70000i
  };
   
  multi_array A(extents[NEL][NMO]);
  multi_array B(extents[NMO][NEL]);

  for(int i=0, k=0; i<A.shape()[0]; i++)
    for(int j=0; j<A.shape()[1]; j++,k++)
       A[i][j] = m_a[k];

  for(int i=0, k=0; i<B.shape()[0]; i++)
    for(int j=0; j<B.shape()[1]; j++,k++)
       B[i][j] = m_b[k];

  multi_array_ref Aref(m_a.data(),extents[NEL][NMO]);
  multi_array_ref Bref(m_b.data(),extents[NMO][NEL]);

  SlaterDetOperations<Type> SDet(NMO,NEL);

  /**** Overlaps ****/
  myREQUIRE(SDet.Overlap(A,B,node),ov);
  myREQUIRE(SDet.Overlap(Aref,B,node),ov);
  myREQUIRE(SDet.Overlap(A,Bref,node),ov);
  myREQUIRE(SDet.Overlap(Aref,Bref,node),ov);

  // Test array_view
  myREQUIRE(SDet.Overlap(A[indices[range_t()][range_t()]],B,node),ov);
  myREQUIRE(SDet.Overlap(A,B[indices[range_t()][range_t()]],node),ov);

  multi_array A_ = A[indices[range_t(0,2)][range_t(0,3)]];
  multi_array B_ = B[indices[range_t(0,3)][range_t(0,2)]];
  myREQUIRE(SDet.Overlap(A[indices[range_t(0,2)][range_t(0,3)]],
                       B[indices[range_t(0,3)][range_t(0,2)]],node),ov2);
  myREQUIRE(SDet.Overlap(A[indices[range_t(0,2)][range_t(0,3)]],B_),ov2);
  myREQUIRE(SDet.Overlap(A_,B[indices[range_t(0,3)][range_t(0,2)]],node),ov2);

  shared_communicator node_ = node.split(node.rank()%2); 
  myREQUIRE(SDet.Overlap(A,B,node_),ov);
  myREQUIRE(SDet.Overlap(A[indices[range_t(0,2)][range_t(0,3)]],
                       B[indices[range_t(0,3)][range_t(0,2)]],node_),ov2);

  /**** Density Matrices *****/
  vector v_ref = {
   1.17573619385025996 - 0.01580426445014660i,  -0.25295981756593167 + 0.28594469607401085i,
  -0.07724502823533341 - 0.09687959052155870i,   0.30512858581808422 - 0.04506898328729603i,

   0.17912592806889663 + 0.08374315906672802i,   0.59381451118048767 + 0.13438888951771200i,
  -0.02021475320201610 - 0.13737561982083193i,   0.32095003919745313 + 0.12832750636154097i,

  -0.04919549564646425 + 0.05402065222741825i,  -0.00286990878355775 - 0.15806420733175885i,
   1.05069081188635494 + 0.00793429988912356i,  -0.08048239150997794 + 0.09917405634760490i,

  -0.46434598219794548 - 0.09896890422706731i,   0.87746748807670427 - 0.53417787485950319i,
   0.12162735438647077 + 0.31042401735800573i,   0.17975848308289613 - 0.12651892495668898i
  };
  vector vc_ref = {
  -0.8721879971495297 + 0.6787593377585239i,   0.3244278932250768 - 1.8083537898275881i,
   0.6130713546530860 + 0.0399736955598931i,   0.0132562806444336 - 0.2882495766584950i,

   0.8853626557603183 + 0.0978868569224204i,  -0.0598704127345155 - 0.0470889603064014i,
  -0.7392693424168300 - 0.0715317395994149i,   0.3721269544963505 - 0.1797896522886788i,

  -0.0567190307022984 - 0.3114847157576828i,  -0.1290126440468128 + 0.6815705660308808i,
   0.3787012855335005 + 0.0039188686237135i,  -0.2005543456941538 + 0.2100886953142371i
  };
  vector v_ref_2 = {
   0.7361983496013835 - 0.0956505507662245i,   0.6467449689807925 + 0.2297471806893873i,
   0.0189270005620390 - 0.1727975708935829i,

   0.2986893588843604 - 0.0099955730815030i,   0.2696068479051557 + 0.0292039710860386i,
   0.0497734835066391 + 0.1783200397050796i,

   0.1020030246826092 + 0.0344707468383766i,  -0.2500340021988402 - 0.0826863644855427i,
   0.9941948024934623 + 0.0664465796801866i
  };
  vector vc_ref_2 = {
  -0.489369975192701 + 0.103038673040713i,  -0.858850485405126 - 0.275734238124941i,
   1.219791948842170 + 0.118626922447301i,
   0.486337033653898 - 0.098631569047656i,   0.595497821653010 + 0.185288949671560i,
  -0.455373025165330 - 0.129360996228044i
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

  SDet.MixedDensityMatrix(A[indices[range_t(0,2)][range_t(0,3)]],
                          B[indices[range_t(0,3)][range_t(0,2)]],
                          G[indices[range_t(0,3)][range_t(0,3)]],node,false);
  check(G[indices[range_t(0,3)][range_t(0,3)]],g_ref_2);
  SDet.MixedDensityMatrix(A_,
                          B[indices[range_t(0,3)][range_t(0,2)]],
                          G[indices[range_t(0,3)][range_t(0,3)]],node,false);
  check(G[indices[range_t(0,3)][range_t(0,3)]],g_ref_2);
  SDet.MixedDensityMatrix(A[indices[range_t(0,2)][range_t(0,3)]],
                          B_,
                          G[indices[range_t(0,3)][range_t(0,3)]],node,false);
  check(G[indices[range_t(0,3)][range_t(0,3)]],g_ref_2);

  SDet.MixedDensityMatrix(A,B,Gc,node,true); check(Gc,gc_ref);
  SDet.MixedDensityMatrix(Aref,B,Gc,node,true); check(Gc,gc_ref);
  SDet.MixedDensityMatrix(A,Bref,Gc,node,true); check(Gc,gc_ref);
  SDet.MixedDensityMatrix(Aref,Bref,Gc,node,true); check(Gc,gc_ref);

  SDet.MixedDensityMatrix(A[indices[range_t(0,2)][range_t(0,3)]],
                          B[indices[range_t(0,3)][range_t(0,2)]],
                          Gc[indices[range_t(0,2)][range_t(0,3)]],node,true);
  check(Gc[indices[range_t(0,2)][range_t(0,3)]],gc_ref_2);
  SDet.MixedDensityMatrix(A_,
                          B[indices[range_t(0,3)][range_t(0,2)]],
                          Gc[indices[range_t(0,2)][range_t(0,3)]],node,true);
  check(Gc[indices[range_t(0,2)][range_t(0,3)]],gc_ref_2);
  SDet.MixedDensityMatrix(A[indices[range_t(0,2)][range_t(0,3)]],
                          B_,
                          Gc[indices[range_t(0,2)][range_t(0,3)]],node,true);
  check(Gc[indices[range_t(0,2)][range_t(0,3)]],gc_ref_2);

  SHM_Buffer SMbuff2(node_,NMO*(NMO+NEL));

  multi_array_ref G2(SMbuff2.data(),extents[NMO][NMO]);
  multi_array_ref Gc2(SMbuff2.data()+NMO*NMO,extents[NEL][NMO]);

  // switch comm
  SDet.MixedDensityMatrix(A,B,G2,node_,false); check(G2,g_ref);
  SDet.MixedDensityMatrix(A[indices[range_t(0,2)][range_t(0,3)]],
                          B[indices[range_t(0,3)][range_t(0,2)]],
                          G2[indices[range_t(0,3)][range_t(0,3)]],node_,false);
  check(G2[indices[range_t(0,3)][range_t(0,3)]],g_ref_2);
  SDet.MixedDensityMatrix(A,B,Gc2,node_,true); check(Gc2,gc_ref);
  SDet.MixedDensityMatrix(A[indices[range_t(0,2)][range_t(0,3)]],
                          B[indices[range_t(0,3)][range_t(0,2)]],
                          Gc2[indices[range_t(0,2)][range_t(0,3)]],node_,true);
  check(Gc2[indices[range_t(0,2)][range_t(0,3)]],gc_ref_2);  

}

TEST_CASE("SDetOps_complex_csr", "[sdet_ops]")
{

  Communicate *c;
  OHMMS::Controller->initialize(0, NULL);
  
  using boost::mpi3::shared_communicator;

  shared_communicator node = boost::mpi3::world.split_shared(boost::mpi3::world.rank()); 

  const int NMO = 4;
  const int NEL = 3;

  using Type = std::complex<double>;
  using vector = std::vector<Type>;
  using multi_array = boost::multi_array<Type,2>;
  using multi_array_ref = boost::multi_array_ref<Type,2>;
  using namespace std::complex_literals;
  using csr_matrix = ma::sparse::csr_matrix<Type,int,int,
                                boost::mpi3::intranode::allocator<Type>,
                                boost::mpi3::intranode::is_root>;

  const Type ov = -7.62332599999999 + 22.20453200000000i;
  const Type ov2 = -10.37150000000000 -  7.15750000000000i;

  // some arbitrary matrices
  vector m_a = {
   0.90000 + 0.10000i,   2.40000 + 0.20000i,   3.00000 + 0.30000i,
   0.40000 + 0.40000i,   1.00000 + 0.50000i,   1.20000 + 0.10000i,
   1.40000 + 0.20000i,   1.60000 + 0.30000i,   3.60000 + 0.40000i,
   0.40000 + 0.50000i,   0.20000 + 0.10000i,   0.10000 + 0.20000i
  };
  vector m_b = {
   1.90000 + 0.60000i,   1.40000 + 0.70000i,   0.40000 + 0.80000i,
   1.40000 + 0.90000i,   0.20000 + 0.50000i,   2.20000 + 0.60000i,
   0.40000 + 0.70000i,   2.60000 + 0.80000i,   0.60000 + 0.90000i,
   1.10000 + 0.50000i,   0.30000 + 0.60000i,   0.90000 + 0.70000i
  };
   
  multi_array A(extents[NMO][NEL]); // Will be transposed when Acsr is built
  multi_array B(extents[NMO][NEL]);

  for(int i=0, k=0; i<A.shape()[0]; i++)
    for(int j=0; j<A.shape()[1]; j++,k++)
       A[i][j] = m_a[k];

  for(int i=0, k=0; i<B.shape()[0]; i++)
    for(int j=0; j<B.shape()[1]; j++,k++)
       B[i][j] = m_b[k];

  multi_array_ref Bref(m_b.data(),extents[NMO][NEL]);

  csr_matrix Acsr(csr::shm::construct_csr_matrix_single_input<csr_matrix>(A,0.0,'T',node));

  SlaterDetOperations<Type> SDet(NMO,NEL);

  /**** Overlaps ****/
  myREQUIRE(SDet.Overlap(Acsr,B,node),ov);
  myREQUIRE(SDet.Overlap(Acsr,Bref,node),ov);

  myREQUIRE(SDet.Overlap(Acsr,B),ov);
  myREQUIRE(SDet.Overlap(Acsr,Bref),ov);

  // Test array_view
  myREQUIRE(SDet.Overlap(Acsr,B[indices[range_t()][range_t()]],node),ov);
  myREQUIRE(SDet.Overlap(Acsr,B[indices[range_t()][range_t()]]),ov);

  shared_communicator node_ = node.split(node.rank()%2); 
  myREQUIRE(SDet.Overlap(Acsr,B,node_),ov);

  multi_array B_ = B[indices[range_t(0,3)][range_t(0,2)]];

  myREQUIRE(SDet.Overlap(Acsr[{0,2,0,3}],B_),ov2);

  /**** Density Matrices *****/
  vector v_ref = {
   1.17573619385025996 - 0.01580426445014660i,  -0.25295981756593167 + 0.28594469607401085i,
  -0.07724502823533341 - 0.09687959052155870i,   0.30512858581808422 - 0.04506898328729603i,

   0.17912592806889663 + 0.08374315906672802i,   0.59381451118048767 + 0.13438888951771200i,
  -0.02021475320201610 - 0.13737561982083193i,   0.32095003919745313 + 0.12832750636154097i,

  -0.04919549564646425 + 0.05402065222741825i,  -0.00286990878355775 - 0.15806420733175885i,
   1.05069081188635494 + 0.00793429988912356i,  -0.08048239150997794 + 0.09917405634760490i,

  -0.46434598219794548 - 0.09896890422706731i,   0.87746748807670427 - 0.53417787485950319i,
   0.12162735438647077 + 0.31042401735800573i,   0.17975848308289613 - 0.12651892495668898i
  };
  vector vc_ref = {
  -0.8721879971495297 + 0.6787593377585239i,   0.3244278932250768 - 1.8083537898275881i,
   0.6130713546530860 + 0.0399736955598931i,   0.0132562806444336 - 0.2882495766584950i,

   0.8853626557603183 + 0.0978868569224204i,  -0.0598704127345155 - 0.0470889603064014i,
  -0.7392693424168300 - 0.0715317395994149i,   0.3721269544963505 - 0.1797896522886788i,

  -0.0567190307022984 - 0.3114847157576828i,  -0.1290126440468128 + 0.6815705660308808i,
   0.3787012855335005 + 0.0039188686237135i,  -0.2005543456941538 + 0.2100886953142371i
  };
  vector v_ref_2 = {
   0.7361983496013835 - 0.0956505507662245i,   0.6467449689807925 + 0.2297471806893873i,
   0.0189270005620390 - 0.1727975708935829i,

   0.2986893588843604 - 0.0099955730815030i,   0.2696068479051557 + 0.0292039710860386i,
   0.0497734835066391 + 0.1783200397050796i,

   0.1020030246826092 + 0.0344707468383766i,  -0.2500340021988402 - 0.0826863644855427i,
   0.9941948024934623 + 0.0664465796801866i
  };
  vector vc_ref_2 = {
  -0.489369975192701 + 0.103038673040713i,  -0.858850485405126 - 0.275734238124941i,
   1.219791948842170 + 0.118626922447301i,
   0.486337033653898 - 0.098631569047656i,   0.595497821653010 + 0.185288949671560i,
  -0.455373025165330 - 0.129360996228044i
  };

  multi_array_ref g_ref(v_ref.data(),extents[NMO][NMO]);
  multi_array_ref gc_ref(vc_ref.data(),extents[NEL][NMO]);
  multi_array_ref g_ref_2(v_ref_2.data(),extents[3][3]);
  multi_array_ref gc_ref_2(vc_ref_2.data(),extents[2][3]);

  using SHM_Buffer = mpi3_SHMBuffer<Type>;
  SHM_Buffer SMbuff(node,NMO*(NMO+NEL));

  multi_array_ref G(SMbuff.data(),extents[NMO][NMO]);
  multi_array_ref Gc(SMbuff.data()+NMO*NMO,extents[NEL][NMO]);

  SDet.MixedDensityMatrix(Acsr,B,G,node,false); check(G,g_ref);
  SDet.MixedDensityMatrix(Acsr,Bref,G,node,false); check(G,g_ref);

  SDet.MixedDensityMatrix(Acsr[{0,2,0,3}],
                          B[indices[range_t(0,3)][range_t(0,2)]],
                          G[indices[range_t(0,3)][range_t(0,3)]],node,false);
  check(G[indices[range_t(0,3)][range_t(0,3)]],g_ref_2);
  SDet.MixedDensityMatrix(Acsr[{0,2,0,3}],
                          B_,
                          G[indices[range_t(0,3)][range_t(0,3)]],node,false);
  check(G[indices[range_t(0,3)][range_t(0,3)]],g_ref_2);

  SDet.MixedDensityMatrix(Acsr,B,Gc,node,true); check(Gc,gc_ref);
  SDet.MixedDensityMatrix(Acsr,Bref,Gc,node,true); check(Gc,gc_ref);

  SDet.MixedDensityMatrix(Acsr[{0,2,0,3}],
                          B[indices[range_t(0,3)][range_t(0,2)]],
                          Gc[indices[range_t(0,2)][range_t(0,3)]],node,true);
  check(Gc[indices[range_t(0,2)][range_t(0,3)]],gc_ref_2);
  SDet.MixedDensityMatrix(Acsr[{0,2,0,3}],
                          B_,
                          Gc[indices[range_t(0,2)][range_t(0,3)]],node,true);
  check(Gc[indices[range_t(0,2)][range_t(0,3)]],gc_ref_2);

  SHM_Buffer SMbuff2(node_,NMO*(NMO+NEL));

  multi_array_ref G2(SMbuff2.data(),extents[NMO][NMO]);
  multi_array_ref Gc2(SMbuff2.data()+NMO*NMO,extents[NEL][NMO]);

  // switch comm
  SDet.MixedDensityMatrix(Acsr,B,G2,node_,false); check(G2,g_ref);
  SDet.MixedDensityMatrix(Acsr[{0,2,0,3}],
                          B[indices[range_t(0,3)][range_t(0,2)]],
                          G2[indices[range_t(0,3)][range_t(0,3)]],node_,false);
  check(G2[indices[range_t(0,3)][range_t(0,3)]],g_ref_2);
  SDet.MixedDensityMatrix(Acsr,B,Gc2,node_,true); check(Gc2,gc_ref);
  SDet.MixedDensityMatrix(Acsr[{0,2,0,3}],
                          B[indices[range_t(0,3)][range_t(0,2)]],
                          Gc2[indices[range_t(0,2)][range_t(0,3)]],node_,true);
  check(Gc2[indices[range_t(0,2)][range_t(0,3)]],gc_ref_2);

}

}
