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

/*! \fn CrystalLattice::CrystalLattice()
 *  Default constructor. Initialized to a \p 1x1x1 cubic supercell.
 */
template<class T, unsigned D>
CrystalLattice<T,D>::CrystalLattice() {
  R.diagonal(1e10);  
  G = R;
  M = R;
  Volume = 1;
  for(int i=0; i<D; i++) BoxBConds[i] = 1;
  reset();
}


template<class T, unsigned D>
CrystalLattice<T,D>::CrystalLattice(const CrystalLattice<T,D>& rhs) {

  R = rhs.R;
  BConds = rhs.BConds;
  for(int idir=0; idir<D; idir++) {
    BoxBConds[idir] = rhs.BoxBConds[idir];//!< Bonundary Condition
  }
  reset();
}

template<class T, unsigned D>
void CrystalLattice<T,D>::set(int argc, char **argv) {
  vector<string> opt;
  for(int i=0; i<argc; i++) opt.push_back(argv[i]);
  set(opt);
}

template<class T, unsigned D>
void CrystalLattice<T,D>::set(vector<string>& argv) {
  makelattice<CrystalLattice<T,D> >::apply(*this, argv);
}



template<class T, unsigned D>
void CrystalLattice<T,D>::set(const Tensor<T,D>& lat) {
   R = lat;
   reset();
}

template<class T, unsigned D>
void CrystalLattice<T,D>::set(T sc, T* lat) {
  if(lat) {
    int ii=0;   
    for(int i=0; i<D; i++)
      for(int j=0; j<D; j++)
        R(i,j) = lat[ii++];
 
    R *= sc;
  } else {
    R = 0.0e0;
    R.diagonal(sc);
  }
  reset();
}

template<class T, unsigned D>
void 
CrystalLattice<T,D>::set(const CrystalLattice<T,D>& oldlat, int *uc) {
  R = oldlat.R;
  BConds = oldlat.BConds;
  for(int idir=0; idir<D; idir++) {
   BoxBConds[idir] = oldlat.BoxBConds[idir];//!< Bonundary Condition
  }
  if(uc) {
    for(int i=0; i<D; i++)
      for(int j=0; j<D; j++) R(i,j) *= static_cast<T>(uc[i]);
  }
  reset();
}

template<class T, unsigned D>
void CrystalLattice<T,D>::reset() {

  G = inverse(R); //G = transpose(Inverse(R));
  Volume = fabs(det(R));
  //M = dot(transpose(R),R);
  M = dot(R,transpose(R));
  T t = pow(TWOPI,2.0);
  Mg = t*dot(transpose(G),G);
  for(int i=0; i<D; i++)
    for(int j=0; j<D; j++)
      Rv[i][j] = R(i,j);
  for(int i=0; i<D; i++)
    for(int j=0; j<D; j++)
      Gv[i][j] = G(j,i);
}

/*! \fn  CrystalLattice<T,D>::operator=(const CrystalLattice<T,D>& rhs)
 *  \param rhs a CrystalLattice to be copied
 *  \brief Copy all the properties of the lattice, 
 *   including boundary conditions and grid paritions.
 */
template<class T, unsigned D>
CrystalLattice<T,D>& 
CrystalLattice<T,D>::operator=(const CrystalLattice<T,D>& rhs) {

  R = rhs.R;
  BConds = rhs.BConds;
  for(int idir=0; idir<D; idir++) {
   BoxBConds[idir] = rhs.BoxBConds[idir];//!< Bonundary Condition
  }
  reset();
  return *this;
}

/*! \fn  CrystalLattice<T,D>::operator=(const Tensor<T,D>& rhs)
 *  \param rhs a Tensor to be copied
 */
template<class T, unsigned D>
CrystalLattice<T,D>& 
CrystalLattice<T,D>::operator=(const Tensor<T,D>& rhs) {
  R = rhs;
  reset();
  return *this;
}

/*! \fn  CrystalLattice<T,D>::operator*=(T sc)
 *  \param sc A scaling factor.
 *  \brief Rescale this supercell by a scalar.
 */
template<class T, unsigned D>
CrystalLattice<T,D>& 
CrystalLattice<T,D>::operator*=(T sc) {
  R *= sc;
  reset();
  return *this;
}

template<class T, unsigned D>
void CrystalLattice<T,D>::print(ostream& os, int level) const {

  /*\note level == 0: print only the lattice vectors
   *      level == 1: lattice vectors, boundary conditions, grid 
   *      level == 2: + all the internal values
   */

  os << "<lattice>" << endl; 
  for(int i=0; i<D; i++) {
    os << Rv[i] << endl;
  }
  os << "</lattice>" << endl; 

  if(level > 0) {
    os << "<bconds> ";
    for(int i=0; i<D; i++) {
      if(BoxBConds[i]) os << " p "; 
      else             os << " n ";
    }
    os << " </bconds>" << endl;
  }

  if(level > 1) {
    os << "Volume (A^3) = " << Volume << endl;
    os << "Reciprocal vectors without 2*pi.\n";
    for(int i=0; i<D; i++) {
      os << "g_"<< i+1<< " = " << Gv[i] << endl;
    }
    os << "Metric tensor in real-space.\n";
    for(int i=0; i<D; i++) {
      os << "h_"<< i+1<< " = ";
      for(int j=0; j< D; j++) {
        os << M(i,j) << " ";
      }
      os << endl;
    }
    os << "Metric tensor in g-space.\n";
    for(int i=0; i<D; i++) {
      os << "h_"<< i+1<< " = ";
      for(int j=0; j< D; j++) {
        os << Mg(i,j) << " ";
      }
      os << endl;
    }
  }
}


template<class T, unsigned D>
inline bool operator==(const CrystalLattice<T,D>& lhs, 
                       const CrystalLattice<T,D>& rhs) {
  const double eps = 1e-6;
  for(int i=0; i<D*D; i++) 
    if(abs(lhs.R[i]-rhs.R[i]) > eps) return false;
  return true;
}

template<class T, unsigned D>
inline bool operator!=(const CrystalLattice<T,D>& lhs, 
                       const CrystalLattice<T,D>& rhs) {
  return !(lhs == rhs);
}


// free function to check if a CrystalLattice is orthorhombic
template<class T, unsigned D>
inline bool orthorombic(const CrystalLattice<T,D>& a) {
  return true;
}

//!< Returns true if the off-diagonal elements of Rm are zero.
template<class T>
bool orthorombic(const CrystalLattice<T,2>& a, T) {
  const T eps = 1e-6;
  return ((a.R(0,1) < eps) && (a.R(1,0) < eps));
}

//!< Returns true if the off-diagonal elements of Rm are zero.
template<class T>
bool orthorombic(const CrystalLattice<T,3>& a, T) {
  const T eps = 1e-6;
  return ((a.R(0,1) < eps) && (a.R(0,2) < eps) &&
          (a.R(1,0) < eps) && (a.R(1,2) < eps) &&
          (a.R(2,0) < eps) && (a.R(2,1) < eps) );
}


/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

