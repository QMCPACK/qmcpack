//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: 
//
// File created by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//////////////////////////////////////////////////////////////////////////////////////
// -*- C++ -*-
/** @file diff_ylm.cpp
 * @brief Check correctness of SoaSphericalTensor and CartesianTensor
 */
#include <Configuration.h>
#include <random/random.hpp>
#include <QMCWaveFunctions/lcao/SoaSphericalTensor.h>
#include <QMCWaveFunctions/lcao/SoaCartesianTensor.h>
#include <Numerics/SphericalTensor.h> 
#include <Numerics/CartesianTensor.h> 
#include <getopt.h>

using namespace std;
using namespace qmcplusplus;

int main(int argc, char** argv)
{

  OHMMS::Controller->initialize(argc,argv);
  if (OHMMS::Controller->rank() != 0) {
    outputManager.shutOff();
  }

  Communicate* mycomm=OHMMS::Controller;

  typedef float RealType;
  typedef TinyVector<RealType,3> PosType;

  //typedef QMCTraits::RealType           RealType;
  //typedef ParticleSet::ParticlePos_t    ParticlePos_t;
  //typedef ParticleSet::ParticleLayout_t LatticeType;
  //typedef ParticleSet::TensorType       TensorType;
  //typedef ParticleSet::PosType          PosType;
  //use the global generator

  bool ionode=(mycomm->rank() == 0);
  int na=4;
  int lmax=4; 
  int nsamples=5;
  char *g_opt_arg;
  int opt;
  while((opt = getopt(argc, argv, "hl:s:")) != -1)
  {
    switch(opt)
    {
      case 'h':
        printf("[-l lmax -s samples]\n");
        return 1;
      case 'l':
        lmax=atoi(optarg);
        break;
      case 's':
        nsamples=atoi(optarg);
        break;
    }
  }

  RandomGenerator<RealType> random(MakeSeed(0,1));
  constexpr RealType small=std::numeric_limits<RealType>::epsilon();

  constexpr RealType shift(0.5);
  TinyVector<double,4> err;

  int ntot;

  //Test SphericalTensor
  {
    SphericalTensor<RealType,PosType> st_0(lmax);
    SoaSphericalTensor<RealType> st_1(lmax);
    ntot=st_0.size();

    for(int i=0; i<nsamples; ++i)
    {
      RealType x=random()-shift;
      RealType y=random()-shift;
      RealType z=random()-shift;
      RealType d=std::sqrt(x*x+y*y+z*z);
      PosType p(x/d,y/d,z/d);

      //test evaluate value
      st_0.evaluate(p);
      st_1.evaluateV(p[0],p[1],p[2]);

      const RealType* ylm=st_1[0];
      for(int lm=0; lm<ntot; ++lm)
      {
        err[0]+= std::abs(st_0.Ylm[lm]-ylm[lm]);
      }

      st_0.evaluateAll(p);
      st_1.evaluateVGL(p[0],p[1],p[2]);

      const RealType* ylm_x=st_1[1];
      const RealType* ylm_y=st_1[2];
      const RealType* ylm_z=st_1[3];

      for(int lm=0; lm<ntot; ++lm)
      {
        PosType ylm_grad(ylm_x[lm],ylm_y[lm],ylm_z[lm]);
        PosType delta_grad=ylm_grad-st_0.gradYlm[lm];
        err[0]+= std::abs(st_0.Ylm[lm]-ylm[lm]);
        err[1]+= std::sqrt(dot(delta_grad,delta_grad));
      }
    }
  }


  //Test cartesian tensor
  if(lmax>4)
   cout << "Skip Cartesian tensor tenor as lmax>4" << endl;
  else
  {
    CartesianTensor<RealType,PosType> cart_0(lmax);
    SoaCartesianTensor<RealType> cart_1(lmax);

    for(int i=0; i<nsamples; ++i)
    {
      RealType x=random()-shift;
      RealType y=random()-shift;
      RealType z=random()-shift;
      RealType d=std::sqrt(x*x+y*y+z*z);
      PosType p(x/d,y/d,z/d);

      cart_0.evaluate(p);
      cart_1.evaluateV(p[0],p[1],p[2]);

      const RealType* xyz=cart_1[0];

      for(int lm=0; lm<ntot; ++lm)
      {
        err[2]+= std::abs(cart_0.XYZ[lm]-xyz[lm]);
      }

      cart_0.evaluateAll(p);
      cart_1.evaluateVGL(p[0],p[1],p[2]);

      const RealType* xyz_x=cart_1[1];
      const RealType* xyz_y=cart_1[2];
      const RealType* xyz_z=cart_1[3];
      for(int lm=0; lm<ntot; ++lm)
      {
        PosType xyz_grad(xyz_x[lm],xyz_y[lm],xyz_z[lm]);
        PosType delta_grad=xyz_grad-cart_0.gradXYZ[lm];
        err[2]+= std::abs(cart_0.XYZ[lm]-xyz[lm]);
        err[3]+= std::sqrt(dot(delta_grad,delta_grad));
      }
    }
  }

  cout << "Error " << ntot << " " << nsamples << " " <<err << endl;

  OHMMS::Controller->finalize();

  return 0;
}
