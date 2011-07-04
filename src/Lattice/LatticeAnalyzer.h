//////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Jeongnim Kim and Kenneth P. Esler
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
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
#ifndef QMCPLUSPLUS_LATTICE_ANALYZER_H
#define QMCPLUSPLUS_LATTICE_ANALYZER_H
#include "OhmmsPETE/TinyVector.h"
namespace qmcplusplus
{

  /** enumeration to clssify a CrystalLattice
  */
  enum {SUPERCELL_OPEN=0, SUPERCELL_WIRE=1, SUPERCELL_SLAB=3, SUPERCELL_BULK=7};

  ///generic class to analyze a Lattice
  template<typename T,  unsigned D> struct LatticeAnalyzer { };

  /** specialization for  3D lattice
  */
  template<typename T>
    struct LatticeAnalyzer<T,3>
    { 

      typedef TinyVector<T,3> SingleParticlePos_t;
      typedef Tensor<T,3> Tensor_t;
      ///SuperCell type
      int mySC;

      inline int operator()(const TinyVector<int,3>& box) 
      {
        return mySC=box[0]+2*(box[1]+box[2]*2);
      }

      inline bool isDiagonalOnly(const Tensor_t& R) const
      {
        T offdiag=abs(R(0,1))+abs(R(0,2))+abs(R(1,0))+abs(R(1,2))+abs(R(2,0))+abs(R(2,1));
        return (offdiag< numeric_limits<T>::epsilon());
      }

      inline SingleParticlePos_t calcSolidAngles(
          const TinyVector<SingleParticlePos_t,3>& Rv, 
          const SingleParticlePos_t& OneOverLength)
      {
        const T rad_to_deg = 180.0/M_PI;
        return SingleParticlePos_t(
            rad_to_deg*std::acos(dot(Rv[0],Rv[1])*OneOverLength[0]*OneOverLength[1]),
            rad_to_deg*std::acos(dot(Rv[1],Rv[2])*OneOverLength[1]*OneOverLength[2]),
            rad_to_deg*std::acos(dot(Rv[2],Rv[0])*OneOverLength[2]*OneOverLength[0]));
      }

      inline T calcWignerSeitzRadius(TinyVector<SingleParticlePos_t,3>& a)
      {
        T rMin = 0.5*std::numeric_limits<T>::max();
        if(mySC == SUPERCELL_BULK) //bulk type
        {
          for (int i=-1; i<=1; i++)
            for (int j=-1; j<=1; j++)
              for (int k=-1; k<=1; k++) 
                if (i || j || k)
                {
                  SingleParticlePos_t L = (static_cast<T>(i) * a[0] + static_cast<T>(j) * a[1] + static_cast<T>(k) * a[2]);
                  rMin=std::min(rMin,dot(L,L));
                }
        }
        else if(mySC == SUPERCELL_SLAB)//slab type
        {
          for (int i=-1; i<=1; i++)
            for (int j=-1; j<=1; j++)
              if (i || j)
              {
                SingleParticlePos_t L = (static_cast<T>(i) * a[0] + static_cast<T>(j) * a[1]);
                rMin=std::min(rMin,dot(L,L));
              }
          cout << " calcWignerSeitzRadius for Slab" << endl;
        }
        else if(mySC == SUPERCELL_WIRE)//wire
        {
          rMin=dot(a[0],a[0]);
        }
        return 0.5*std::sqrt(rMin);
      }

      inline T calcSimulationCellRadius(TinyVector<SingleParticlePos_t,3>& a)
      {
        T scr = 0.5*std::numeric_limits<T>::max();
        //if(mySC == SUPERCELL_BULK)
        //{
          for (int i=0; i<3; ++i) 
          {
            SingleParticlePos_t A = a(i);
            SingleParticlePos_t B = a((i+1)%3);    
            SingleParticlePos_t C = a((i+2)%3);
            SingleParticlePos_t BxC = cross(B,C);
            T dist = 0.5*std::abs(dot(A,BxC))/std::sqrt(dot(BxC,BxC));
            scr = std::min(scr, dist);
          }
        //}
        //else if(mySC == SUPERCELL_SLAB)
        //{
        //  T a0mag  = std::sqrt(dot(a[0],a[0]));
        //  T a1mag  = std::sqrt(dot(a[1],a[1]));
        //  scr=0.5*std::min(a0mag,a1mag);
        //  //T dist = 0.5*std::abs(dot(A,BxC))/std::sqrt(dot(BxC,BxC));
        //  //T theta1 = dot (a[0], a[1])/(a0mag*a1mag);
        //  //T theta2 = M_PI - theta1;
        //  //T theta  = std::min (theta1, theta2);
        //  //T dist   = std::min (a0mag, a1mag);
        //  //scr=0.5*std::sin(theta)*dist;
        //  cout << " calcSimulationCellRadius for Slab" << endl;
        //}
        //else if(mySC == SUPERCELL_WIRE)
        //{
        //  scr=0.5*std::sqrt(dot(a[0],a[0]));
        //}
        return scr;
      }

      inline void makeNextCells(const Tensor_t& lat, vector<SingleParticlePos_t>& nextcells)
      {
        int ic=0;
        if(mySC == SUPERCELL_BULK)
        {
          nextcells.resize(26);
          for(int i=-1;i<=1;++i)
            for(int j=-1;j<=1;++j)
              for(int k=-1;k<=1;++k)
              {
                if(!(i || j || k )) continue;//exclude zero
                SingleParticlePos_t u(i,j,k);
                nextcells[ic++]=DotProduct<SingleParticlePos_t,Tensor_t,false>::apply(u,lat);    
              }
        }
        else if(mySC == SUPERCELL_SLAB)
        {
          nextcells.resize(8);
          int k=0;
          for(int i=-1;i<=1;++i)
            for(int j=-1;j<=1;++j)
              if(i||j)
              {
                SingleParticlePos_t u(i,j,k);
                nextcells[ic++]=DotProduct<SingleParticlePos_t,Tensor_t,false>::apply(u,lat);    
              }
        }
        else if(mySC == SUPERCELL_WIRE)
        {
          nextcells.resize(2);
          SingleParticlePos_t um(-1,0,0);
          nextcells[ic++]=DotProduct<SingleParticlePos_t,Tensor_t,false>::apply(um,lat);    
          SingleParticlePos_t up(1,0,0);
          nextcells[ic++]=DotProduct<SingleParticlePos_t,Tensor_t,false>::apply(up,lat);    
        }
      }
    };

  /** specialization for  2D lattice
  */
  template<typename T>
    struct LatticeAnalyzer<T,2>
    { 
      typedef TinyVector<T,2> SingleParticlePos_t;
      typedef Tensor<T,2> Tensor_t;

      /** return supercell enum
       * @param[in] box[2] if box[i]==1, PBC
       * @return SUPERCELL_OPEN or SUPERCELL_BULK
       */
      inline int operator()(const TinyVector<int,2>& box) 
      {
        return (box[0]+2*box[1])? SUPERCELL_BULK: SUPERCELL_OPEN;
      }

      inline bool isDiagonalOnly(const Tensor<T,2>& R) const
      {
        T offdiag=abs(R(0,1))+abs(R(1,0));
        return (offdiag< numeric_limits<T>::epsilon());
      }

      inline SingleParticlePos_t calcSolidAngles(
          const TinyVector<SingleParticlePos_t,2>& Rv, 
          const SingleParticlePos_t& OneOverLength)
      {
        const T rad_to_deg = 180.0/M_PI;
        return SingleParticlePos_t(
            rad_to_deg*std::acos(dot(Rv[0],Rv[1])*OneOverLength[0]*OneOverLength[1])
            , 0.0);
      }

      inline T calcWignerSeitzRadius(TinyVector<SingleParticlePos_t,2>& a)
      {
        T rMin = 1.0e50;
        for (int i=-1; i<=1; i++)
          for (int j=-1; j<=1; j++)
            if ((i!=0) || (j!=0)) {
              SingleParticlePos_t L = ((double)i * a[0] + (double)j * a[1]);
              T dist = 0.5*std::fabs(dot(L,L));
              rMin = std::min(rMin, dist);
            }
        return rMin;
      }

      inline T calcSimulationCellRadius(TinyVector<SingleParticlePos_t,2>& a)
      {
	T a0mag  = std::sqrt(dot(a[0],a[0]));
	T a1mag  = std::sqrt(dot(a[1],a[1]));
	T theta1 = dot (a[0], a[1])/(a0mag*a1mag);
        T theta2 = M_PI - theta1;
	T theta  = std::min (theta1, theta2);
	T dist   = std::min (a0mag, a1mag);
	return 0.5*std::sin(theta)*dist;
        // return calcWignerSeitzRadius(a);
      }
      
      inline void makeNextCells(const Tensor_t& lat, vector<SingleParticlePos_t>& nextcells)
      {
        nextcells.resize(8);
        int ic=0;
        for(int i=-1;i<=1;++i)
          for(int j=-1;j<=1;++j)
            {
              if(!(i || j)) continue;//exclude zero
              SingleParticlePos_t u(i,j);
              nextcells[ic++]=DotProduct<SingleParticlePos_t,Tensor_t,false>::apply(u,lat);    
            }
      }
      
    };

  /** specialization for 1D lattice
  */
  template<typename T>
    struct LatticeAnalyzer<T,1>
    { 
      typedef TinyVector<T,1> SingleParticlePos_t;
      inline bool isDiagonalOnly(const Tensor<T,1>& R) const
      {
        return true;
      }

      inline int operator()(const TinyVector<int,1>& box) 
      {
        return (box[0])? SUPERCELL_BULK:SUPERCELL_OPEN;
      }

      inline T calcWignerSeitzRadius(TinyVector<SingleParticlePos_t,1>& a)
      {
        return a[0]*0.5;
      }
    };

}
#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1189 $   $Date: 2006-07-17 10:08:11 -0500 (Mon, 17 Jul 2006) $
 * $Id: DistanceTable.cpp 1189 2006-07-17 15:08:11Z jnkim $ 
 ***************************************************************************/
