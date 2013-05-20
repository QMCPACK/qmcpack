//////////////////////////////////////////////////////////////////
// (c) Copyright 1998-2002 by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//   Department of Physics, Ohio State University
//   Ohio Supercomputer Center
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
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
  std::vector<string> opt;
  for(int i=0; i<argc; i++)
    opt.push_back(argv[i]);
  set(opt);
}

template<class T, unsigned D,bool ORTHO>
void CrystalLattice<T,D,ORTHO>::set(vector<string>& argv)
{
  makelattice<CrystalLattice<T,D,ORTHO> >::apply(*this, argv);
}


template<class T, unsigned D,bool ORTHO>
void CrystalLattice<T,D,ORTHO>::set(const Tensor<T,D>& lat)
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
  //analysis of a lattice using LatticeAnalyzer
  LatticeAnalyzer<T,D> ldesc;
  SuperCellEnum = ldesc(BoxBConds);
  DiagonalOnly=ldesc.isDiagonalOnly(R);
  ABC=ldesc.calcSolidAngles(Rv,OneOverLength);
  WignerSeitzRadius = ldesc.calcWignerSeitzRadius(Rv);
  SimulationCellRadius = ldesc.calcSimulationCellRadius(Rv);
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
  //  T dist = 0.5*std::fabs(dot(A,BxC))/std::sqrt(dot(BxC,BxC));
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
void CrystalLattice<T,D,ORTHO>::print(ostream& os, int level) const
{
  /*\note level == 0: print only the lattice vectors
   *      level == 1: lattice vectors, boundary conditions, grid
   *      level == 2: + all the internal values
   */
  os << "<parameter name=\"lattice\">" << endl;
  for(int i=0; i<D; ++i)
    os << Rv[i] << endl;
  os << "</parameter>" << endl;
  if(level > 0)
  {
    os << "<parameter name=\"bconds\"> ";
    for(int i=0; i<D; ++i)
    {
      if(BoxBConds[i])
        os << " p ";
      else
        os << " n ";
    }
    os << "</parameter>" << endl;
  }
  os << "<note>"<<endl;
  if(level > 1)
  {
    os << "Volume (A^3) = " << Volume << endl;
    os << "Reciprocal vectors without 2*pi.\n";
    for(int i=0; i<D; ++i)
      os << "g_"<< i+1<< " = " << Gv[i] << endl;
    os << "Metric tensor in real-space.\n";
    for(int i=0; i<D; ++i)
    {
      os << "h_"<< i+1<< " = ";
      for(int j=0; j< D; ++j)
      {
        os << M(i,j) << " ";
      }
      os << endl;
    }
    os << "Metric tensor in g-space.\n";
    for(int i=0; i<D; ++i)
    {
      os << "h_"<< i+1<< " = ";
      for(int j=0; j< D; ++j)
      {
        os << Mg(i,j) << " ";
      }
      os << endl;
    }
  }
  os << "</note>"<<endl;
}

template<class T, unsigned D,bool ORTHO>
inline bool operator==(const CrystalLattice<T,D,ORTHO>& lhs,
                       const CrystalLattice<T,D,ORTHO>& rhs)
{
  for(int i=0; i<D*D; ++i)
    if(abs(lhs.R[i]-rhs.R[i]) > numeric_limits<T>::epsilon())
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

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/

