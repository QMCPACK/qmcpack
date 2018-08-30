//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2018 QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
/** \file SplineAdoptorVectorized.h
 * \class SplineAdoptorVectorized
 * Base class for SplineAdoptor's used for BsplineSet<SplineAdoptor>
 * Specifies that a SplineXXXAdopter provides these functions
 * - evaluate_v    value only
 * - evaluate_vgl  vgl
 * - evaluate_vgh  vgh
 * Specializations are implemented  in Spline*Adoptor.h and include
 * - SplineC2RAdoptor<ST,TT,D> : real wavefunction using complex einspline, tiling
 * - SplineC2CAdoptor<ST,TT,D> : complex wavefunction using complex einspline, tiling
 * - SplineR2RAdoptor<ST,TT,D> : real wavefunction using real einspline, a single twist
 * where ST (TT) is the precision of the einspline (SPOSetBase).
 *
 * typedefs and data members are duplicated for each adoptor class.
 * @todo Specalization and optimization for orthorhombic cells to use vgl not vgh
 */
#ifndef QMCPLUSPLUS_SPLINEADOPTORBASEVECTORIZED_H
#define QMCPLUSPLUS_SPLINEADOPTORBASEVECTORIZED_H

//#include "QMCWaveFunctions/BsplineFactory/BsplineDevice.h"
//#include "QMCWaveFunctions/BsplineFactory/BsplineDeviceCUDA.h"
#include "Lattice/CrystalLattice.h"
#include "simd/allocator.hpp"

namespace qmcplusplus
{

/** base class any SplineAdoptor
 *
 * This handles SC and twist and declare storage for einspline
 */
template<template<typename, unsigned> class DEVICE, typename ST, unsigned D>
class SplineAdoptorVectorized
{
public:
  //static_assert(std::is_base_of<BsplineDevice<DEVICE, ST, D>, DEVICE>, "DEVICE must inherit from BsplineDevice");
  using PointType = TinyVector<ST,D>;
  using SingleSplineType = typename DEVICE<ST,D>::SingleBsplineType;
  using DataType = ST; 

  DEVICE<ST,D> bspline_dev;
  ///true if the computed values are complex
  bool is_complex;
  ///true, if it has only one k point and Gamma
  bool is_gamma_only;
  ///true, if it has only one k point and Gamma
  bool is_soa_ready;
  ///Index of this adoptor, when multiple adoptors are used for NUMA or distributed cases
  size_t MyIndex;
  ///number of unique orbitals
  size_t nunique_orbitals;
  ///first index of the SPOs this Spline handles
  size_t first_spo;
  ///last index of the SPOs this Spline handles
  size_t last_spo;
  ///sign bits at the G/2 boundaries
  TinyVector<int,D>          HalfG;
  ///\f$GGt=G^t G \f$, transformation for tensor in LatticeUnit to CartesianUnit, e.g. Hessian
  Tensor<ST,D>               GGt;
  CrystalLattice<ST,D>       SuperLattice;
  CrystalLattice<ST,D>       PrimLattice;
  /// flags to unpack sin/cos
  std::vector<bool>               MakeTwoCopies;
  ///kpoints for each unique orbitals
  std::vector<TinyVector<ST,D> >  kPoints;
  ///remap band
  aligned_vector<int> BandIndexMap;
  /// band offsets
  std::vector<int> offset_real, offset_cplx;
  ///name of the adoptor
  std::string AdoptorName;
  ///keyword used to match hdf5
  std::string KeyWord;
  SplineAdoptorVectorized()
    :is_complex(false),is_gamma_only(false), is_soa_ready(false),
    MyIndex(0),nunique_orbitals(0),first_spo(0),last_spo(0)
  { }
  SplineAdoptorVectorized(const SplineAdoptorVectorized& rhs)=default;

  inline void init_base(int n)
  {
    nunique_orbitals=n;
    GGt=dot(transpose(PrimLattice.G),PrimLattice.G);
    kPoints.resize(n);
    MakeTwoCopies.resize(n);
  }

  ///remap kpoints to group general kpoints & special kpoints
  int remap_kpoints()
  {
    std::vector<TinyVector<ST,D> >  k_copy(kPoints);
    const int nk=kPoints.size();
    BandIndexMap.resize(nk);
    int nCB=0;
    //two pass
    for(int i=0; i<nk; ++i)
    {
      if(MakeTwoCopies[i]) 
      {
        kPoints[nCB]=k_copy[i];
        BandIndexMap[i]=nCB++;
      }
    }
    int nRealBands=nCB;
    for(int i=0; i<nk; ++i)
    {
      if(!MakeTwoCopies[i]) 
      {
        kPoints[nRealBands]=k_copy[i];
        BandIndexMap[i]=nRealBands++;
      }
    }
    return nCB; //return the number of complex bands
  }
};

}
#endif
