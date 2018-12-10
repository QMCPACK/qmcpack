//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



/**@file CrystalLattice.cpp
 *@brief Member function definitions of CrystalLattice<T,D>, included by CrystalLattice.h
 */
#include "Lattice/MakeCrystalLattice.h"
#include "Lattice/LatticeAnalyzer.h"

namespace qmcplusplus
{
/*! \fn CrystalLattice::CrystalLattice()
 *  Default constructor. Initialized to a \p 1x1x1 cubic supercell.
 */
template<class T, unsigned D,bool ORTHO>
CrystalLattice<T,D,ORTHO>::CrystalLattice()
{
  BoxBConds=0;
  VacuumScale=1.0;
  R.diagonal(1e10);
  G = R;
  M = R;
  Volume = 1;
  reset();
}

//  template<class T, unsigned D,bool ORTHO>
//    CrystalLattice<T,D,ORTHO>::CrystalLattice(const CrystalLattice<T,D>& rhs)
//    {
//      BoxBConds = rhs.BoxBConds;
//      R = rhs.R;
//      reset();
//    }

template<class T, unsigned D,bool ORTHO>
void CrystalLattice<T,D,ORTHO>::set(int argc, char **argv)
{
  std::vector<std::string> opt;
  for(int i=0; i<argc; i++)
    opt.push_back(argv[i]);
  set(opt);
}

template<class T, unsigned D,bool ORTHO>
void CrystalLattice<T,D,ORTHO>::set(std::vector<std::string>& argv)
{
  makelattice<CrystalLattice<T,D,ORTHO> >::apply(*this, argv);
}


template<class T, unsigned D,bool ORTHO>
template<class TT>
void CrystalLattice<T,D,ORTHO>::set(const Tensor<TT,D>& lat)
{
  R = lat;
  reset();
}

template<class T, unsigned D,bool ORTHO>
void CrystalLattice<T,D,ORTHO>::set(T sc, T* lat)
{
  if(lat)
  {
    for(int i=0; i<D; ++i)
      for(int j=0; j<D; ++j)
        R(i,j) = *lat++;
    R *= sc;
  }
  else
  {
    R = 0.0e0;
    R.diagonal(sc);
  }
  reset();
}

template<class T, unsigned D,bool ORTHO>
void
CrystalLattice<T,D,ORTHO>::set(const CrystalLattice<T,D,ORTHO>& oldlat, int *uc)
{
  BoxBConds = oldlat.BoxBConds;
  VacuumScale = oldlat.VacuumScale;
  R = oldlat.R;
  if(uc)
  {
    for(int i=0; i<D; ++i)
      for(int j=0; j<D; ++j)
        R(i,j) *= static_cast<T>(uc[i]);
  }
  reset();
}

template<class T, unsigned D,bool ORTHO>
void CrystalLattice<T,D,ORTHO>::reset()
{
  G = inverse(R); //G = transpose(Inverse(R));
  Gt= transpose(G);
  Volume = std::abs(det(R));
  //M = dot(transpose(R),R);
  M = dot(R,transpose(R));
  T t = TWOPI*TWOPI;
  Mg = t*dot(transpose(G),G);
  for(int i=0; i<D; ++i)
    for(int j=0; j<D; ++j)
      Rv[i][j] = R(i,j);
  for(int i=0; i<D; ++i)
    for(int j=0; j<D; ++j)
      Gv[i][j] = G(j,i);
  for(int i=0; i<D; ++i)
  {
    Length[i]=std::sqrt(dot(Rv[i],Rv[i]));
    OneOverLength[i]=1.0/Length[i];
  }
  Center = 0.0;
  for(int i=0;i<D;++i)
    Center += Rv[i];
  Center *= .5;
  //analysis of a lattice using LatticeAnalyzer
  LatticeAnalyzer<T,D> ldesc;
  SuperCellEnum = ldesc(BoxBConds);
  DiagonalOnly=ldesc.isDiagonalOnly(R);
  ABC=ldesc.calcSolidAngles(Rv,OneOverLength);
  WignerSeitzRadius = ldesc.calcWignerSeitzRadius(Rv);
  WignerSeitzRadius_G = ldesc.calcWignerSeitzRadius(Gv);
  SimulationCellRadius = ldesc.calcSimulationCellRadius(Rv);
  // set equal WignerSeitzRadius and SimulationCellRadius when they are very close.
  if ( WignerSeitzRadius > SimulationCellRadius &&
       WignerSeitzRadius-SimulationCellRadius <= WignerSeitzRadius*std::numeric_limits<float>::epsilon()*2 )
    WignerSeitzRadius = SimulationCellRadius;
  CellRadiusSq=SimulationCellRadius*SimulationCellRadius;
  if(SuperCellEnum)
    ldesc.makeNextCells(R,NextUnitCells);
  //SimulationCellRadius = 1.0e50;
  //// Compute simulation cell radius
  //for (int i=0; i<D; i++) {
  //  SingleParticlePos_t A = a(i);
  //  SingleParticlePos_t B = a((i+1)%3);
  //  SingleParticlePos_t C = a((i+2)%3);
  //  SingleParticlePos_t BxC = cross(B,C);
  //  T dist = 0.5*std::abs(dot(A,BxC))/std::sqrt(dot(BxC,BxC));
  //  SimulationCellRadius = std::min(SimulationCellRadius, dist);
  //}
}

//  /*! \fn  CrystalLattice<T,D>::operator=(const CrystalLattice<T,D>& rhs)
//   *  \param rhs a CrystalLattice to be copied
//   *  \brief Copy all the properties of the lattice,
//   *   including boundary conditions and grid paritions.
//   */
//  template<class T, unsigned D,bool ORTHO>
//    CrystalLattice<T,D,ORTHO>&
//    CrystalLattice<T,D,ORTHO>::operator=(const CrystalLattice<T,D,ORTHO>& rhs)
//    {
//      BoxBConds = rhs.BoxBConds;
//      R = rhs.R;
//      reset();
//      return *this;
//    }
//
///*! \fn  CrystalLattice<T,D>::operator=(const Tensor<T,D>& rhs)
// *  \param rhs a Tensor to be copied
// */
//template<class T, unsigned D,bool ORTHO>
//  CrystalLattice<T,D,ORTHO>&
//  CrystalLattice<T,D,ORTHO>::operator=(const Tensor<T,D>& rhs)
//  {
//    R = rhs;
//    reset();
//    return *this;
//  }

/*! \fn  CrystalLattice<T,D>::operator*=(T sc)
 *  \param sc A scaling factor.
 *  \brief Rescale this supercell by a scalar.
 */
template<class T, unsigned D,bool ORTHO>
CrystalLattice<T,D,ORTHO>&
CrystalLattice<T,D,ORTHO>::operator*=(T sc)
{
  R *= sc;
  reset();
  return *this;
}

template<class T, unsigned D,bool ORTHO>
void CrystalLattice<T,D,ORTHO>::print(std::ostream& os, int level) const
{
  /*\note level == 0: print only the lattice vectors
   *      level == 1: lattice vectors, boundary conditions, grid
   *      level == 2: + all the internal values
   */
  std::string unit_name = "bohr";

  std::string lattice_name = "  Lattice (" + unit_name + "):";
  std::string pad(lattice_name.length(),' ');
  os <<  lattice_name;
  for(int i=0; i<D; ++i) {
    if (i > 0) {
      os << pad;
    }
    os << Rv[i] << std::endl;
  }
  if(level > 0)
  {
    os << std::endl;
    os << "  Boundary Conditions: ";
    for(int i=0; i<D; ++i)
    {
      if(BoxBConds[i])
        os << " p ";
      else
        os << " n ";
    }
    os << std::endl;
    if(VacuumScale != 1.0)
      os << "  Vacuum scale: " << VacuumScale << std::endl;
  }
  if(level > 1)
  {
    os << std::endl;
    os << "  Volume (bohr^3) = " << Volume << std::endl;
    os << std::endl;
    os << "  Reciprocal vectors without 2*pi.\n";
    for(int i=0; i<D; ++i)
      os << "    g_"<< i+1<< " = " << Gv[i] << std::endl;
    os << std::endl;
    os << "  Metric tensor in real-space.\n";
    for(int i=0; i<D; ++i)
    {
      os << "    h_"<< i+1<< " = ";
      for(int j=0; j< D; ++j)
      {
        os << M(i,j) << " ";
      }
      os << std::endl;
    }
    os << std::endl;
    os << "  Metric tensor in g-space.\n";
    for(int i=0; i<D; ++i)
    {
      os << "    h_"<< i+1<< " = ";
      for(int j=0; j< D; ++j)
      {
        os << Mg(i,j) << " ";
      }
      os << std::endl;
    }
  }
}

template<class T, unsigned D,bool ORTHO>
inline bool operator==(const CrystalLattice<T,D,ORTHO>& lhs,
                       const CrystalLattice<T,D,ORTHO>& rhs)
{
  for(int i=0; i<D*D; ++i)
    if(std::abs(lhs.R[i]-rhs.R[i]) > std::numeric_limits<T>::epsilon())
      return false;
  return true;
}

template<class T, unsigned D,bool ORTHO>
inline bool operator!=(const CrystalLattice<T,D,ORTHO>& lhs,
                       const CrystalLattice<T,D,ORTHO>& rhs)
{
  return !(lhs == rhs);
}

// free function to check if a CrystalLattice is orthorhombic
template<class T, unsigned D,bool ORTHO>
inline bool orthorombic(const CrystalLattice<T,D,ORTHO>& a)
{
  return ORTHO;
}
}


