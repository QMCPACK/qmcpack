//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#ifndef GUARD_HETEROSTRUCTURE_H
#define GUARD_HETEROSTRUCTURE_H

#include "Numerics/Spline3D/Config.h"
#include "Numerics/Spline3D/Grid3D.h"
#include "Numerics/Spline3D/Spline3D.h"
#include "Particle/ParticleBase.h"
#include "QMCHamiltonians/WOS/Layer.h"
#include "QMCHamiltonians/WOS/Interface.h"
#include "QMCHamiltonians/WOS/Domain.h"


class HeteroStructure
{

public:
  int nlayers_d;

  /// skin width
  double skin_d;

  /// weight due to interfaces
  double weight_bc;

  /// different factors
  double Gfac;

  /// boundaries of HeteroStructure
  posvec_t r_min;
  posvec_t r_max;

  /// HeteroStructure on a Grid3D
  Grid3D* Full_Grid;

  /// Spline for the charge density
  Spline3D* rho_Spline;

  std::vector<double> eps_d;
  std::vector<double> offset_d;
  std::vector<double> intervals_d;

  std::vector<Layer*> layers;
  std::vector<Interface*> interfaces;

  std::vector<double> bare_potential;
  std::vector<double> rhoq;




  /// initialise and construct the HeteroStructure
  HeteroStructure(double delta,
                  Grid3D* aGrid3D,
                  const char* rhofile,
                  const char* vfile,
                  std::vector<double>& dielectric,
                  std::vector<double>& BandOffset,
                  std::vector<double>& allIntervals)
  {
    // Gfac = e/4pieps0*nm = -1.6e-19 * 9e9 / 1e-9
    Gfac = -1.439952241;
    skin_d = delta;
    Full_Grid = aGrid3D;
    rho_Spline = new Spline3D(aGrid3D,aGrid3D->npts_m);
    /// number of layers is one less than intervals
    nlayers_d = allIntervals.size() - 1;
    for( int i = 0; i < nlayers_d; i++)
    {
      eps_d.push_back(dielectric[i]);
      offset_d.push_back(BandOffset[i]);
      intervals_d.push_back(allIntervals[i]);
    }
    intervals_d.push_back(allIntervals[nlayers_d]);
    /// initialise r_min and r_max
    for( int dir = 0; dir < 3; dir++ )
    {
      r_min[ dir ] = Full_Grid->d_axis[ dir ].d_x[ 0 ];
      r_max[ dir ] = Full_Grid->d_axis[ dir ].d_x[ Full_Grid->nvec_m[dir]-1];
    }
    /// initilaise charge and potential vectors
    bare_potential.resize(Full_Grid->size);
    rhoq.resize(Full_Grid->size);
    /// construct layers and interfaces
    construct();
    /// read in the files
    read_file( -1, vfile, bare_potential );
    read_file( 1, rhofile, rhoq );
    /// read in the charge into the Spline3D
    //rho_Spline->read_data(rhofile);
    //rho_Spline->d2fdr2();
    //    for(int i = 0; i < rhoq.size(); i++)
    //      bare_potential[i] *= 0.036749033500418936;
  }

  /// initialise the weight for the applied potential
  inline void weight_init()
  {
    weight_bc = 1.0;
  }

  /// construct layers and interfaces
  void construct();

  /// Spherical Domain contained within layer
  void MakeLimitedSphere( Domain& );

  /// Spherical Domain centered on an interface
  void MakeIfcSphere( Domain&, int );

  /// Maximum Spherical Domain entirely contained within HeteroStructure
  void MakeMaximumSphere( Domain& );

  /// add a new layer to HeteroStructure
  void add_Layer( Layer* );

  /// add a new interface to HeteroStructure
  void add_Interface( Interface* );

  /// calculate dfrac
  void calc_dfrac( Domain& );

  /// probability of sampling each section of domain
  void sample_prob( const Domain& );

  /// dielectric at position r
  double epsilon( const posvec_t& );

  /// dielectric at domain center :: runner
  void epsilon( Domain& );

  /// sample a point on the Domian surface
  void sample_point( Domain& );

  /// collect contribution from charge density + other point charges Gee
  double Image_Contribution( int,
                             const Domain&, const qmcplusplus::ParticleBase& );

  /// collect contribution from charge density + other point charges Gee
  double Domain_Contribution( const Domain&, const qmcplusplus::ParticleBase& );

  /// read in the boundary potential and charge density
  void read_file( int, const char*, std::vector<double>& );

  /// first passage
  double passage( Domain& );



};

#endif
