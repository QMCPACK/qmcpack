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
    
    



#ifndef GUARD_TRICUBICSPLINESET_H
#define GUARD_TRICUBICSPLINESET_H

#include <vector>
#include <map>
#include <libxml/xmlmemory.h>
#include <libxml/parser.h>

#include "Numerics/Spline3D/Config.h"
#include "Numerics/Spline3D/Grid3D.h"
#include "Numerics/Spline3D/SetSplinePoint.h"
#include "Numerics/Spline3D/TriCubicSpline.h"

namespace ohmmsqmc
{

class TriCubicSplineSet
{

  bool OwnGrid;

  /// the number of Spline orbitals
  int norbs;

  /// the initial and final corners of the Schroedinger's region
  gridvec_t ri,rf;

  /// pointer to the full grid
  Grid3D* DeviceGrid;

  /// pointer to the Schroedinger grid
  Grid3D* Schr_Grid;

  /// the vector of pointers to all the orbitals
  std::vector<TriCubicSpline*> m_psi;

  /// class to set r-point parameters
  SetSplinePoint* m_set;

  /// read all the wave-functions from the input file
  void readwf(const char*, int ispin, int num);

  std::map<std::string,int> orbital_map;

public:

  /// constructor
  TriCubicSplineSet();

  /// destructor
  ~TriCubicSplineSet();

  /// reading in the input data
  bool put(xmlNodePtr, Grid3D*);


  /// a function to set the two grid pointers
  inline void set(Grid3D* fg, Grid3D* Sg)
  {
    DeviceGrid = fg;
    Schr_Grid = Sg;
  }

  inline int size() const
  {
    return m_psi.size();
  }

  /// return the fullGrid
  Grid3D* get_DeviceGrid()
  {
    return DeviceGrid;
  }

  /// return an orbital
  TriCubicSpline* getOrbital(int i)
  {
    return m_psi[i];
  }

  /// setting the point r parameters
  inline void set_point(const posvec_t& r)
  {
    m_set->set_point(r,Schr_Grid);
  }

  /// initialising the padded function
  void finit(int);

  inline TriCubicSpline* getOrbital(const std::string& orbname)
  {
    std::map<std::string,int>::iterator it = orbital_map.find(orbname);
    if(it == orbital_map.end())
      return NULL;
    else
      return m_psi[(*it).second];
  }

  /// evaluate the wavefunction
  void evaluate(const posvec_t&,scalar_array_t&,posarray_t&,scalar_array_t&);

  /// more esoteric evaluation
  template<class PTCL, class VM, class GM>
  inline void evaluate(const PTCL& P, int first, int last,
                       VM& logdet, GM& dlogdet, VM& d2logdet)
  {
    int n = last-first;
    int iat = first;
    for(int i=0; i < n; i++,iat++)
    {
      set_point(P.R[iat],Schr_Grid);   /// set the r-point parameters
      for(int j=0; j < n; j++)
      {
        logdet(j,i)= m_psi[j]->evaluate(P.R[iat], dlogdet(i,j),d2logdet(i,j));
      }
    }
  }





};

}
#endif
