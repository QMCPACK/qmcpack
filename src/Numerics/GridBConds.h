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
    
    



/** @file GridBConds.h
 * @brief Definitions of GridBCond classes
 */
#ifndef QMCPLUSPLUS_GRIDBC_HANDLER_H
#define QMCPLUSPLUS_GRIDBC_HANDLER_H

#include <config/stdlib/math.h>

namespace qmcplusplus
{

/** enumeration for the template parameter BC to handle GridBCHanlder
 */
enum {FIXED_GBC=0, PERIODIC_GBC, NO_GBC};

/** general grid handler
 *
 * boolean member is_periodic is used to swith on/off periodic boundary conditions
 */
template<typename T, unsigned BC=NO_GBC>
struct GridBCond
{
  ///true, if periodic
  bool periodic;
  /// number of grid points
  int Ng;
  /// lower bound of a grid
  T Lower;
  /// upper bound of a grid
  T Upper;
  /// 1/Delta
  T OneOverDelta;
  /// Upper-Lower
  T L;
  /// 1.0/(Upper-Lower)
  T OneOverL;
  ///grid spacing,  (Upper-Lower)/Ng
  T Delta;

  ///default constructor
  inline GridBCond(): periodic(true),Lower(0.0), Upper(1.0) {}

  ///return periodic
  inline bool is_periodic() const
  {
    return periodic;
  }

  /** initialize the grid
   * @param ng number of grid points
   * @param xmin lower bound
   * @param xmax upper bound
   */
  inline void init(int ng, T xmin, T xmax, bool pbc)
  {
    periodic=pbc;
    Ng=ng;
    Lower=xmin;
    Upper=xmax;
    L=xmax-xmin;
    OneOverL=1.0/L;
    Delta=L/static_cast<T>(Ng);
    OneOverDelta=1.0/Delta;
  }

  /** apply boundary condition to a difference
   * @param x a separation of two points
   */
  inline void applyBC(T& x) const
  {
    if(periodic)
    {
      T x1=x*OneOverL;
      x=L*(x1-round(x1));
    }
  }

  /** evaluate the grid index for a position x
   * @param x position, shifted by Lower
   * @param i index of the shifted position i
   * @return false, if x = [Lower,Upper)
   */
  inline bool outofbound(T& x, int& i) const
  {
    x -=Lower;
    if(periodic)
      x -= std::floor(x*OneOverL)*L;
    else
      if(x<0.0 || x>=L)
      {
        return true;
      }
    T xi;
    x = modf (x*OneOverDelta, &xi);
    i = static_cast<int>(xi);
    return false;
  }
};

/** handler of an one-dimensional grid with vanising boundary conditions
 *
 * One dimensional grid is defined on [Lower,Upper).
 */
template<typename T>
struct GridBCond<T,FIXED_GBC>
{
  /// number of grid points
  int Ng;
  /// lower bound of a grid
  T Lower;
  /// upper bound of a grid
  T Upper;
  /// 1/Delta
  T OneOverDelta;
  /// Upper-Lower
  T L;
  /// 1.0/(Upper-Lower)
  T OneOverL;
  ///grid spacing,  (Upper-Lower)/Ng
  T Delta;

  ///default constructor
  inline GridBCond(): Lower(0.0), Upper(1.0) {}

  ///return false
  inline bool is_periodic() const
  {
    return false;
  }

  /** initialize the grid
   * @param ng number of grid points
   * @param xmin lower bound
   * @param xmax upper bound
   * @param pbc ignored
   */
  inline void init(int ng, T xmin, T xmax, bool pbc)
  {
    Ng=ng;
    Lower=xmin;
    Upper=xmax;
    L=xmax-xmin;
    OneOverL=1.0/L;
    Delta=L/static_cast<T>(Ng);
    OneOverDelta=1.0/Delta;
  }

  /** apply boundary condition to a difference
   * @param x a separation of two points
   */
  inline void applyBC(T& x) const {}

  /** evaluate the grid index for a position x
   * @param x position, shifted by Lower
   * @param i index of the shifted position i
   * @return false, if x = [Lower,Upper)
   */
  inline bool outofbound(T& x, int& i) const
  {
    if(x<Lower || x>=Upper)
    {
      return true;
    }
    T xi;
    x = modf ((x-Lower)*OneOverDelta, &xi);
    i = static_cast<int>(xi);
    return false;
  }
};

/** handler of an one-dimensional grid with periodic boundary conditions */
template<typename T>
struct GridBCond<T,PERIODIC_GBC>
{
  int Ng;
  T Lower;
  T L;
  T OneOverL;
  T OneOverDelta;
  T Delta;
  T Upper;

  inline GridBCond(): L(1.0),OneOverL(1.0) {}

  ///return true
  inline bool is_periodic() const
  {
    return true;
  }

  inline void init(int ng, T xmin, T xmax, bool pbc)
  {
    Ng=ng;
    Lower=xmin;
    Upper=xmax;
    L=xmax-xmin;
    OneOverL=1.0/L;
    Delta=L/static_cast<T>(Ng);
    OneOverDelta=1.0/Delta;
  }

  /** evaluate the grid index for a position x
   * @param x position, shifted by Lower
   * @param i index of the shifted position i
   * @return false
   */
  inline bool outofbound(T& x, int& i) const
  {
    T xi;
    x -= Lower;
    x -= std::floor(x*OneOverL)*L;
    x = modf (x*OneOverDelta, &xi);
    i = static_cast<int>(xi);
    return false;
  }

  inline void applyBC(T& x) const
  {
    T dmy;
    T x1=std::modf(x*OneOverL,&dmy);
    //T x1=std::fmod(x*OneOverL,1.0);
    x=L*(x1-static_cast<int>(2.0*x1));
  }
};
}
#endif
