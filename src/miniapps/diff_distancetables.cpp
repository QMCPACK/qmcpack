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
/** @file distancetables_soa.cpp
 *
 * Test code for the accuracy of AoS to SoA transformation of distance tables.
 */
#include <Configuration.h>
#include <Particle/ParticleSet.h>
#include <Particle/DistanceTable.h>
#include <Particle/DistanceTableData.h>
#include <OhmmsSoA/VectorSoaContainer.h>
#include <random/random.hpp>
#include <mpi/collectives.h>
#include <miniapps/input.hpp>
#include <miniapps/pseudo.hpp>
#include <Utilities/Timer.h>
#include <miniapps/common.hpp>
#include <getopt.h>

using namespace std;
using namespace qmcplusplus;

int main(int argc, char** argv)
{

  OHMMS::Controller->initialize(0, NULL);
  Communicate* mycomm=OHMMS::Controller;

  //use the global generator

  bool ionode=(mycomm->rank() == 0);
  int na=4;
  int nb=4;
  int nc=1;
  int nsteps=1;
  int iseed=11;
  int nx=48,ny=48,nz=60;
  //thread blocking
  //int ncrews=1; //default is 1
  int tileSize=-1;
  int ncrews=1;

  char *g_opt_arg;
  int opt;
  while((opt = getopt(argc, argv, "hs:g:i:")) != -1)
  {
    switch(opt)
    {
      case 'h':
        printf("[-g \"n0 n1 n2\"]\n");
        return 1;
      case 'g': //tiling1 tiling2 tiling3
        sscanf(optarg,"%d %d %d",&na,&nb,&nc);
        break;
      case 's'://random seed
        iseed=atoi(optarg);
        break;
      case 'i'://number of iterations
        nsteps=atoi(optarg);
        break;
    }
  }

  Random.init(0,1,iseed);

  typedef QMCTraits::RealType           RealType;
  typedef ParticleSet::ParticlePos_t    ParticlePos_t;
  typedef ParticleSet::ParticleLayout_t LatticeType;
  typedef ParticleSet::TensorType       TensorType;
  typedef ParticleSet::PosType          PosType;

  Tensor<int,3> tmat(na,0,0,0,nb,0,0,0,nc);
  double t0=0.0,t1=0.0;
  constexpr  OHMMS_PRECISION scale(1);

  RandomGenerator<RealType> random_th(11);

  ParticleSet ions, els;

  ions.setName("ion0");
  ions.Lattice.BoxBConds=1;
  ions.Lattice.LR_rc=5;
  tile_cell(ions,tmat,scale);
  ions.RSoA=ions.R; //this needs to be handled internally

  const int nions=ions.getTotalNum();
  const int nels=count_electrons(ions);
  const int nels3=3*nels;

  {//create up/down electrons
    els.Lattice.BoxBConds=1;   
    els.Lattice.LR_rc=5;
    els.Lattice.set(ions.Lattice);
    vector<int> ud(2); ud[0]=nels/2; ud[1]=nels-ud[0];
    els.create(ud);
    els.R.InUnit=1;
    random_th.generate_uniform(&els.R[0][0],nels3);
    els.convert2Cart(els.R); // convert to Cartiesian
    els.RSoA=els.R; //this needs to be handled internally
    els.setName("e");
  }

  constexpr RealType eps=numeric_limits<float>::epsilon();

  //copy of ParticleSet for validations
  ParticleSet els_aos(els);
  ParticleSet ions_aos(ions);

  ParticleSet::ParticlePos_t Rcopy(els.R);

  DistanceTableData* d_ee=DistanceTable::add(els,DT_SOA);
  DistanceTableData* d_ie=DistanceTable::add(ions,els,DT_SOA);

  d_ie->setRmax(els.Lattice.WignerSeitzRadius);
  RealType Rsim=els.Lattice.WignerSeitzRadius;
  std::cout << "Setting 1-body cutoff " <<  d_ie->Rmax << std::endl;

  DistanceTableData* d_ee_aos=DistanceTable::add(els_aos,DT_AOS);
  DistanceTableData* d_ie_aos=DistanceTable::add(ions_aos,els_aos,DT_AOS);

  //SoA version does not need update if PbyP
  els.update();
  els_aos.update();

  //SoA check symmetry
  double sym_err=0.0;
  double sym_all_err=0.0;
  double sym_disp_err=0.0;
  double sym_disp_all_err=0.0;
  int nn=0;
  for(int iel=0; iel<nels; ++iel)
    for(int jel=iel+1; jel<nels; ++jel)
    {
      RealType dref=d_ee_aos->r(nn);
      RealType dsym=std::abs(d_ee->Distances[jel][iel]-d_ee->Distances[iel][jel]);
      PosType dr= d_ee->Displacements[jel][iel]+d_ee->Displacements[iel][jel];
      RealType d2=sqrt(dot(dr,dr));
      sym_all_err += dsym;
      sym_disp_all_err+=d2;
      if(dref<Rsim)
      {
        sym_err += dsym;
        sym_disp_err+=d2;
      }
      ++nn;
    }
  cout << "---------------------------------" << endl;
  cout << "AA SoA(upper) - SoA(lower) (ALL) distances     = " << sym_all_err/nn << endl;
  cout << "AA SoA(upper) + SoA(lower) (ALL) displacements = " << sym_disp_all_err/nn << endl;
  cout << "AA SoA(upper) - SoA(lower) distances     = " << sym_err/nn << endl;
  cout << "AA SoA(upper) + SoA(lower) displacements = " << sym_disp_err/nn << endl;
  
  ParticlePos_t delta(nels);
  
  //main particle-by-particle update
  RealType sqrttau=0.2;
  for(int s=0; s<nsteps; ++s)
  {
    random_th.generate_normal(&delta[0][0],nels3);
    for(int iel=0; iel<nels; ++iel)
    {
      els.setActive(iel);    
      els_aos.setActive(iel);    

      PosType dr=sqrttau*delta[iel];
      bool valid_move=els.makeMoveAndCheck(iel,dr); 
      valid_move &= els_aos.makeMoveAndCheck(iel,dr); 

      if(valid_move && Random()<0.5)
      {
        els.rejectMove(iel);
        els_aos.rejectMove(iel);
      }
      else
      {
        els.acceptMove(iel);
        els_aos.acceptMove(iel);
      }
    }
  }

  els.donePbyP();
  els_aos.donePbyP();

  std::cout << "WignerSeitzRadius " <<  d_ie->Rmax << std::endl;
  {
    double r_err=0.0;
    for(int iat=0; iat<nels; ++iat)
    {
      PosType d=els.R[iat]-Rcopy[iat];
      RealType r=sqrt(dot(d,d));
      r_err += r;
    }
    cout << "Done with the sweep. Diffusion |els.R-R0|^2/nels = " << r_err/nels << endl;
  }
  {
    double r_err=0.0;
    for(int iat=0; iat<nels; ++iat)
    {
      PosType d=els.R[iat]-els_aos.R[iat];
      RealType r=sqrt(dot(d,d));
      r_err += r;
      if(r > eps) 
        cerr << "SoA-AoS error " << iat << " "  << els.RSoA[iat] << " " << els_aos.R[iat] << endl;
    }
    cout << "Done with the sweep. |els.R-els_aos.R|^2/nels = " << r_err/nels << endl;
  }

  {
    double dist_err=0.0;
    double dist_all_err=0.0;
    double disp_err=0.0;
    double disp_all_err=0.0;
    int nn=0;
    for(int iel=0; iel<nels; ++iel)
      for(int jel=iel+1; jel<nels; ++jel)
      {
        RealType dref=d_ee_aos->r(nn);
        RealType d= std::abs(d_ee->Distances[jel][iel]-dref);
        PosType dr= (d_ee->Displacements[jel][iel]+d_ee_aos->dr(nn));
        RealType d2=sqrt(dot(dr,dr));
        dist_all_err+=d;
        disp_all_err += d2;
        if(dref<Rsim)
        {
          dist_err += d;
          disp_err += d2;
        }
        ++nn;
      }
    cout << "---------------------------------" << endl;
    cout << "AA SoA-AoS in distance     (ALL) = " << dist_all_err/nn << endl;
    cout << "AA SoA-AoS in displacement (ALL) = " << disp_all_err/nn << endl;
    cout << endl;
    cout << "AA SoA-AoS in distance           = " << dist_err/nn << endl;
    cout << "AA SoA-AoS in displacement       = " << disp_err/nn << endl;
  }

  {
    double dist_err=0.0;
    double dist_all_err=0.0;
    double disp_all_err=0.0;
    double disp_err=0.0;
    int nn=0;
    for(int i=0; i<nions; ++i)
      for(int j=0; j<nels; ++j)
      {
        RealType dref=d_ie_aos->r(nn);
        RealType d= std::abs(d_ie->Distances[j][i]-dref);
        PosType dr= (d_ie->Displacements[j][i]+d_ie_aos->dr(nn));
        RealType d2=sqrt(dot(dr,dr));
        dist_all_err += d;
        disp_all_err += d2;
        if(dref<Rsim)
        {
          dist_err += d;
          disp_err += d2;
        }
        ++nn;
      }
    cout << "---------------------------------" << endl;
    cout << "AB SoA-AoS in distance     (ALL) = " << dist_all_err/nn << endl;
    cout << "AB SoA-AoS in displacement (ALL) = " << disp_all_err/nn << endl;
    cout << endl;
    cout << "AB SoA-AoS in distance           = " << dist_err/nn << endl;
    cout << "AB SoA-AoS in displacement       = " << disp_err/nn << endl;
  }


  {
    double dist_err=0.0;
    double displ_err=0.0;
    int nn=0;
    for(int i=0; i<nions; ++i)
      for(int nj=0, jmax=d_ie->M[i]; nj<jmax; ++nj)
      {
        int j=d_ie->J2(i,nj);
        int ij=i*nels+j;
        PosType diff_dr=d_ie->dr_m2(i,nj)-d_ie_aos->dr(ij);
        RealType diff_r=d_ie->r_m2(i,nj)- d_ie_aos->r(ij);
        dist_err += abs(diff_r);
        displ_err += sqrt(dot(diff_dr,diff_dr));
        nn++;
      }
    cout << "---------------------------------" << endl;
    cout << "AB SoA-AoS compact = " << dist_err/nn << " " << displ_err/nn;
    cout << " Total number of nn = " << nn << " / " << nels*nions << endl;
  }

  OHMMS::Controller->finalize();

  return 0;
}
