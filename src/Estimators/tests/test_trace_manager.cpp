//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2018 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Mark Dewing, mdewin@anl.gov, Argonne National Laboratory
//
// File created by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"



#include "Message/Communicate.h"
#include "OhmmsData/Libxml2Doc.h"
// Hack to convert aborts into test failures.
// Works because all abort calls are in the .h file.
#undef APP_ABORT
#define APP_ABORT(msg) FAIL(msg)
#include "Estimators/TraceManager.h"


#include <stdio.h>
#include <sstream>

namespace qmcplusplus
{


// See QMCHamiltonian.cpp::initialize_traces for the usage sequence

TEST_CASE("TraceManager", "[estimators]")
{
  OHMMS::Controller->initialize(0, NULL);
  Communicate *c = OHMMS::Controller;

  Array<TraceInt, 1>* int_sample;

  TraceManager tm(c);

  tm.put(NULL, true, "test");

  TraceRequest req1;
  req1.contribute_scalar("scalar1", true);
  req1.streaming_scalar("scalar1");

  tm.request.incorporate(req1);

  tm.request.determine_stream_write();

  int_sample = tm.checkout_int<1>("scalar1");
  (*int_sample)(0) = 2;

  tm.update_status();
}

// Moved here from TraceManager.cpp, check_trace_build
TEST_CASE("TraceManager check_trace_build", "[estimators]")
{
  OHMMS::Controller->initialize(0, NULL);

  std::string domain = "domain";
  std::string name   = "name";
  std::string name1   = "name_1";
  std::string name2   = "name_2";
  std::string name3   = "name_3";
  std::string name4   = "name_4";
  int index     = 0;
  int dim       = 3;
  std::string label  = "label";
  std::vector<int> vi;
  std::vector<double> vr;
  std::vector<std::complex<double> > vc;
  // Shape checking code requires size > 0
  TinyVector<int,4> shape1 = {1,0,0,0};
  TinyVector<int,4> shape2 = {1,1,0,0};
  TinyVector<int,4> shape3 = {1,1,1,0};
  TinyVector<int,4> shape4 = {1,1,1,1};
  ParticleSet P;
  P.create(1); // zero-sized particle set not handled well by TraceManager
  TraceSample<int> tsi(domain,name,index,dim,vi);
  TraceSample<double> tsr(domain,name,index,dim,vr);
  TraceSample<std::complex<double> > tsc(domain,name,index,dim,vc);
  TraceSamples<int> tssi;
  TraceSamples<TraceReal> tssr;
  TraceSamples<TraceComp> tssc;
  TraceBuffer<int> tbi;
  TraceBuffer<TraceReal> tbr;

  Array<int,1>* ai1;
  Array<int,2>* ai2;
  Array<int,3>* ai3;
  Array<int,4>* ai4;
  ai1 = tssi.checkout_array<1>(domain,name1,shape1);
  ai1 = tssi.checkout_array<1>(P,name,shape1);
  ai2 = tssi.checkout_array<2>(domain,name2,shape2);
  ai2 = tssi.checkout_array<2>(P,name2,shape2);
  ai3 = tssi.checkout_array<3>(domain,name3,shape3);
  ai3 = tssi.checkout_array<3>(P,name3,shape3);
  ai4 = tssi.checkout_array<4>(domain,name4,shape4);
  ai4 = tssi.checkout_array<4>(P,name4,shape4);

  Array<TraceReal,1>* ar1;
  Array<TraceReal,2>* ar2;
  Array<TraceReal,3>* ar3;
  Array<TraceReal,4>* ar4;
  ar1 = tssr.checkout_array<1>(domain,name1,shape1);
  ar1 = tssr.checkout_array<1>(P,name1,shape1);
  ar2 = tssr.checkout_array<2>(domain,name2,shape2);
  ar2 = tssr.checkout_array<2>(P,name2,shape2);
  ar3 = tssr.checkout_array<3>(domain,name3,shape3);
  ar3 = tssr.checkout_array<3>(P,name3,shape3);
  ar4 = tssr.checkout_array<4>(domain,name4,shape4);
  ar4 = tssr.checkout_array<4>(P,name4,shape4);


  Array<TraceComp,1>* ac1;
  Array<TraceComp,2>* ac2;
  Array<TraceComp,3>* ac3;
  Array<TraceComp,4>* ac4;
  ac1 = tssc.checkout_array<1>(domain,name1,shape1);
  ac1 = tssc.checkout_array<1>(P,name1,shape1);
  ac2 = tssc.checkout_array<2>(domain,name2,shape2);
  ac2 = tssc.checkout_array<2>(P,name2,shape2);
  ac3 = tssc.checkout_array<3>(domain,name3,shape3);
  ac3 = tssc.checkout_array<3>(P,name3,shape3);
  ac4 = tssc.checkout_array<4>(domain,name4,shape4);
  ac4 = tssc.checkout_array<4>(P,name4,shape4);

  TraceManager tm;
  Array<long, 1>* al1;
  Array<long, 2>* al2;
  Array<long, 3>* al3;
  Array<long, 4>* al4;
  al1 = tm.checkout_int<1>(name1);
  al1 = tm.checkout_int<1>(domain,name1,10);
  al1 = tm.checkout_int<1>(name1,P);
  al2 = tm.checkout_int<2>(name2,5,6);
  al2 = tm.checkout_int<2>(domain,name2,10,11);
  al2 = tm.checkout_int<2>(name2,P,11);
  al3 = tm.checkout_int<3>(name3,5,6,7);
  al3 = tm.checkout_int<3>(domain,name3,10,11,12);
  al3 = tm.checkout_int<3>(name3,P,11,12);
  al4 = tm.checkout_int<4>(name4,5,6,7,8);
  al4 = tm.checkout_int<4>(domain,name4,10,11,12,13);
  al4 = tm.checkout_int<4>(name4,P,11,12,13);
  ar1 = tm.checkout_real<1>(name1);
  ar1 = tm.checkout_real<1>(domain,name1,10);
  ar1 = tm.checkout_real<1>(name1,P);
  ar2 = tm.checkout_real<2>(name2,5,6);
  ar2 = tm.checkout_real<2>(domain,name2,10,11);
  ar2 = tm.checkout_real<2>(name2,P,11);
  ar3 = tm.checkout_real<3>(name3,5,6,7);
  ar3 = tm.checkout_real<3>(domain,name3,10,11,12);
  ar3 = tm.checkout_real<3>(name3,P,11,12);
  ar4 = tm.checkout_real<4>(name4,5,6,7,8);
  ar4 = tm.checkout_real<4>(domain,name4,10,11,12,13);
  ar4 = tm.checkout_real<4>(name4,P,11,12,13);

  ac1 = tm.checkout_complex<1>(name1);
  ac1 = tm.checkout_complex<1>(domain,name1,10);
  ac1 = tm.checkout_complex<1>(name1,P);
  ac2 = tm.checkout_complex<2>(name2,5,6);
  ac2 = tm.checkout_complex<2>(domain,name2,10,11);
  ac2 = tm.checkout_complex<2>(name2,P,11);
  ac3 = tm.checkout_complex<3>(name3,5,6,7);
  ac3 = tm.checkout_complex<3>(domain,name3,10,11,12);
  ac3 = tm.checkout_complex<3>(name3,P,11,12);
  ac4 = tm.checkout_complex<4>(name4,5,6,7,8);
  ac4 = tm.checkout_complex<4>(domain,name4,10,11,12,13);
  ac4 = tm.checkout_complex<4>(name4,P,11,12,13);

}

}
