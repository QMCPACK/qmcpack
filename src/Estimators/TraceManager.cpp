//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
    
    



#if !defined(DISABLE_TRACEMANAGER)


#include <Estimators/TraceManager.h>


namespace qmcplusplus
{




double TraceManager::trace_tol = 1e-8;



//#define CHECK_TRACE_BUILD

#ifdef CHECK_TRACE_BUILD

//note: the presence of check_trace_build causes a multiple-definition problem at link-time
void check_trace_build()
{
  std::string domain = "domain";
  std::string name   = "name";
  int index     = 0;
  int dim       = 3;
  std::string label  = "label";
  std::vector<int> vi;
  std::vector<double> vr;
  std::vector<std::complex<double> > vc;
  TinyVector<int,4> shape;
  ParticleSet P;
  TraceSample<int> tsi(domain,name,index,dim,vi);
  TraceSample<double> tsr(domain,name,index,dim,vr);
  TraceSample<std::complex<double> > tsc(domain,name,index,dim,vc);
  TraceSamples<int> tssi;
  TraceSamples<double> tssr;
  TraceSamples<std::complex<double> > tssc;
  TraceBuffer<int> tbi;
  TraceBuffer<double> tbr;
  Array<int,1>* ai1;
  Array<int,2>* ai2;
  Array<int,3>* ai3;
  Array<int,4>* ai4;
  ai1 = tssi.checkout_array<1>(domain,name,shape);
  ai1 = tssi.checkout_array<1>(P,name,shape);
  ai2 = tssi.checkout_array<2>(domain,name,shape);
  ai2 = tssi.checkout_array<2>(P,name,shape);
  ai3 = tssi.checkout_array<3>(domain,name,shape);
  ai3 = tssi.checkout_array<3>(P,name,shape);
  ai4 = tssi.checkout_array<4>(domain,name,shape);
  ai4 = tssi.checkout_array<4>(P,name,shape);
  Array<double,1>* ar1;
  Array<double,2>* ar2;
  Array<double,3>* ar3;
  Array<double,4>* ar4;
  ar1 = tssr.checkout_array<1>(domain,name,shape);
  ar1 = tssr.checkout_array<1>(P,name,shape);
  ar2 = tssr.checkout_array<2>(domain,name,shape);
  ar2 = tssr.checkout_array<2>(P,name,shape);
  ar3 = tssr.checkout_array<3>(domain,name,shape);
  ar3 = tssr.checkout_array<3>(P,name,shape);
  ar4 = tssr.checkout_array<4>(domain,name,shape);
  ar4 = tssr.checkout_array<4>(P,name,shape);
  Array<std::complex<double>,1>* ac1;
  Array<std::complex<double>,2>* ac2;
  Array<std::complex<double>,3>* ac3;
  Array<std::complex<double>,4>* ac4;
  ac1 = tssc.checkout_array<1>(domain,name,shape);
  ac1 = tssc.checkout_array<1>(P,name,shape);
  ac2 = tssc.checkout_array<2>(domain,name,shape);
  ac2 = tssc.checkout_array<2>(P,name,shape);
  ac3 = tssc.checkout_array<3>(domain,name,shape);
  ac3 = tssc.checkout_array<3>(P,name,shape);
  ac4 = tssc.checkout_array<4>(domain,name,shape);
  ac4 = tssc.checkout_array<4>(P,name,shape);
  TraceManager tm;
  ai1 = tm.checkout_int<1>(name);
  ai1 = tm.checkout_int<1>(domain,name,10);
  ai1 = tm.checkout_int<1>(P,name);
  ai2 = tm.checkout_int<2>(name,5,6);
  ai2 = tm.checkout_int<2>(domain,name,10,11);
  ai2 = tm.checkout_int<2>(P,name,11);
  ai3 = tm.checkout_int<3>(name,5,6,7);
  ai3 = tm.checkout_int<3>(domain,name,10,11,12);
  ai3 = tm.checkout_int<3>(P,name,11,12);
  ai4 = tm.checkout_int<4>(name,5,6,7,8);
  ai4 = tm.checkout_int<4>(domain,name,10,11,12,13);
  ai4 = tm.checkout_int<4>(P,name,11,12,13);
  ar1 = tm.checkout_real<1>(name);
  ar1 = tm.checkout_real<1>(domain,name,10);
  ar1 = tm.checkout_real<1>(P,name);
  ar2 = tm.checkout_real<2>(name,5,6);
  ar2 = tm.checkout_real<2>(domain,name,10,11);
  ar2 = tm.checkout_real<2>(P,name,11);
  ar3 = tm.checkout_real<3>(name,5,6,7);
  ar3 = tm.checkout_real<3>(domain,name,10,11,12);
  ar3 = tm.checkout_real<3>(P,name,11,12);
  ar4 = tm.checkout_real<4>(name,5,6,7,8);
  ar4 = tm.checkout_real<4>(domain,name,10,11,12,13);
  ar4 = tm.checkout_real<4>(P,name,11,12,13);
  ac1 = tm.checkout_complex<1>(name);
  ac1 = tm.checkout_complex<1>(domain,name,10);
  ac1 = tm.checkout_complex<1>(P,name);
  ac2 = tm.checkout_complex<2>(name,5,6);
  ac2 = tm.checkout_complex<2>(domain,name,10,11);
  ac2 = tm.checkout_complex<2>(P,name,11);
  ac3 = tm.checkout_complex<3>(name,5,6,7);
  ac3 = tm.checkout_complex<3>(domain,name,10,11,12);
  ac3 = tm.checkout_complex<3>(P,name,11,12);
  ac4 = tm.checkout_complex<4>(name,5,6,7,8);
  ac4 = tm.checkout_complex<4>(domain,name,10,11,12,13);
  ac4 = tm.checkout_complex<4>(P,name,11,12,13);
}
#endif

}


#endif

