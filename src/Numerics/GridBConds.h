//////////////////////////////////////////////////////////////////
// (c) Copyright 2006-  Kenneth Esler  and Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Modified by Jeongnim Kim for qmcpack
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
/** @file GridBConds.h
 * @brief Definitions of GridBCond classes
 */
#ifndef QMCPLUSPLUS_GRIDBC_HANDLER_H
#define QMCPLUSPLUS_GRIDBC_HANDLER_H

namespace qmcplusplus {

  /** enumeration for the template parameter BC to handle GridBCHanlder
   */
  enum {FIXED_GBC=0, PERIODDIC_GBC, NO_GBC};

  /** grid handler to be specialized */
  template<typename T, unsigned BC=NO_GBC> 
    struct GridBCond { 
      bool is_periodic;
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
      inline GridBCond(): is_periodic(true),Lower(0.0), Upper(1.0){}

      /** initialize the grid
       * @param ng number of grid points
       * @param xmin lower bound
       * @param xmax upper bound
       */
      inline void init(int ng, T xmin, T xmax, bool pbc)
      {
        is_periodic=pbc;
        Ng=ng;
        Lower=xmin;
        Upper=xmax;
        L=xmax-xmin; 
        OneOverL=1.0/L;
        Delta=L/static_cast<int>(Ng);
        OneOverDelta=1.0/Delta;
      }

      /** apply boundary condition to a difference
       * @param x a separation of two points
       */
      inline void applyBC(T& x) const {
        if(is_periodic)
        {
          T x1=std::fmod(x*OneOverL,1.0);
          x=L*(x1-static_cast<int>(2.0*x1));
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
        if(is_periodic)
          x -= std::floor(x*OneOverL)*L;
        else
          if(x<Lower || x>=Upper) {return true;}
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
    const bool is_periodic=false;
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
    inline GridBCond(): Lower(0.0), Upper(1.0){}

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
      Delta=L/static_cast<int>(Ng);
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
      x -=Lower;
      if(x<Lower || x>=Upper) {return true;}
      T xi;
      x = modf (x*OneOverDelta, &xi);
      i = static_cast<int>(xi);
      return false;
    }
  };

  /** handler of an one-dimensional grid with periodic boundary conditions */
  template<typename T>
  struct GridBCond<T,PERIODIC_GBC>
  {
    const bool is_periodic=true;
    int Ng;
    T Lower;
    T L;
    T OneOverL;
    T OneOverDelta;
    T Delta;
    T Upper;
    inline GridBCond(): L(1.0),OneOverL(1.0){}
    inline void init(int ng, T xmin, T xmax, bool pbc) 
    {
      Ng=ng;
      Lower=xmin;
      Upper=xmax;
      L=xmax-xmin; 
      OneOverL=1.0/L;
      Delta=L/static_cast<int>(Ng);
      OneOverDelta=1.0/Delta;
    }

    /** evaluate the grid index for a position x
     * @param x position, shifted by Lower
     * @param i index of the shifted position i
     * @return false 
     */
    inline bool outofbound(T& x, int& i) const 
    { 
      x-=Lower;
      x -= std::floor(x*OneOverL)*L;
      T xi;
      x = modf (x*OneOverDelta, &xi);
      i = static_cast<int>(xi);
      return false;
    }

    inline void applyBC(T& x) const 
    {
      T x1=std::fmod(x*OneOverL,1.0);
      x=L*(x1-static_cast<int>(2.0*x1));
    }
  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 2109 $   $Date: 2007-06-26 20:53:50 -0500 (Tue, 26 Jun 2007) $
 * $Id: TricubicBsplineGridBC.h 2109 2007-06-27 01:53:50Z jnkim $
 ***************************************************************************/
