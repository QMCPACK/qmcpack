//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#ifndef QMCPLUSPLUS_TRICUBICSPLINEPOTENTIAL_H
#define QMCPLUSPLUS_TRICUBICSPLINEPOTENTIAL_H
#include <algo.h>
#include "Particle/ParticleSet.h"
#include "Particle/WalkerSetRef.h"
#include "Particle/DistanceTableData.h"
#include "Numerics/Spline3D/Grid3D.h"
#include "Numerics/Spline3D/TriCubicSpline.h"
#include "QMCHamiltonians/QMCHamiltonianBase.h"


namespace qmcplusplus
{


struct TriCubicSplinePotential: public QMCHamiltonianBase
{

  double Efac;

  /// pointer to the main grid (initialised in wavefunction)
  Grid3D* DeviceGrid;

  /// the spline to calculate the potential
  TriCubicSpline* pot_m;

  /// Constructor
  TriCubicSplinePotential(double mbyepsq,
                          Grid3D* agrid,
                          const std::string& fname)
  {
    const double Ha = 27.2113845;   /// Hartree in eV
    const double effHa = mbyepsq * Ha; /// effective Hartree
    DeviceGrid = agrid;
    pot_m = new TriCubicSpline(agrid);
    double ufac = 0.036749033500418936;
    Efac = 1.0/effHa;
    /// Create the spline from the given grid and initialise from the file
    std::cout << "Converting Energy units to effective Ha: 1 eV = "
         << Efac << " Ha*, 1Ha* = " << effHa << "" << std::endl;
    std::cout << "Reading Potential File and initialising ... ";
    pot_m->read_data(fname.c_str(),Efac);
    std::cout << "done! " << std::endl;
  }

  /// Destructor
  ~TriCubicSplinePotential() { }

  /// evaluate the potential
  inline ValueType evaluate(ParticleSet& P)
  {
    ValueType e = 0.0;
    for(int i=0; i<P.getTotalNum(); i++)
    {
      pot_m->set_point(P.R[i]);
      e+=pot_m->evaluate(P.R[i]);
      /// offset :: CHANGE WARNING !!!!!
      if(P.R[i][2] < 9.1854)
        e += 0.33*Efac;
    }
    return e;
  }

  inline Return_t evaluate(ParticleSet& P, std::vector<NonLocalData>& Txy)
  {
    return evaluate(P);
  }

  inline ValueType evaluate(ParticleSet& P, RealType& x)
  {
    return x=evaluate(P);
  }

  void evaluate(WalkerSetRef& P, ValueVectorType& LE)
  {
  }

};

}
#endif
