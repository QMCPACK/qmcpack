//////////////////////////////////////////////////////////////////
// (c) Copyright 2006-  Kenneth Esler and Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Modified by Jeongnim Kim for qmcpack
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
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef TRICUBIC_B_SPLINE_SET_H
#define TRICUBIC_B_SPLINE_SET_H

#include "QMCWaveFunctions/OrbitalTraits.h"
#include "Numerics/TricubicBsplineGrid.h"

namespace qmcplusplus {

  template<typename T>
    class TricubicBspline: public OrbitalTraits<T>
    {
      public:
        typedef typename OrbitalTraits<T>::real_type real_type;
        typedef TricubicBsplineGrid<T> GridType;

        TricubicBspline(){}

        inline void setGrid(const GridType& knots)
        {
          bKnots=knots;
        }

        inline void setGrid(real_type xi, real_type xf, 
            real_type yi, real_type yf, real_type zi, real_type zf, 
            int nx, int ny, int nz, 
            bool interp=true, bool periodic=true,bool openend=true)
        {
          bKnots.setGrid(xi,xf,yi,yf,zi,zf,nx,ny,nz,interp,periodic,openend);
        }

        void Init(const Array<T,3>& data)
        {
          bKnots.Init(data,P);
        }

        inline T evaluate(const TinyVector<real_type,3>& r) 
        {
          bKnots.Find(r[0],r[1],r[2]);
          return bKnots.evaluate(P);
        }

        inline T 
          evaluate(const TinyVector<real_type,3>& r, TinyVector<T,3>& gradf, T& lapf)
          {
            bKnots.FindAll(r[0],r[1],r[2]);
            return bKnots.evaluate(P,gradf,lapf);
          }

      private:
        //Grid
        GridType bKnots;
        // The control points
        Array<T,3> P;
    };

  template<typename T>
    class TricubicBsplineSet: public OrbitalTraits<T>
    {
      public:
        typedef typename OrbitalTraits<T>::real_type real_type;
        typedef typename OrbitalTraits<T>::value_type value_type;
        typedef TricubicBsplineGrid<T> GridType;
        typedef Array<T,3>             StorageType;

        /** default constructure
         *
         * For memory efficiency, reserve DeleteP and P
         * OffSet is set to 1000000. Safe until we can do 1000000 orbitals.
         */
        TricubicBsplineSet():OffSet(1000000) 
        { 
          DeleteP.reserve(1024); 
          P.reserve(1024);
        }

        ~TricubicBsplineSet() { 
          for(int i=0; i<DeleteP.size(); i++)
          {
            if(DeleteP[i]) delete P[i];
          }
        }

        inline void setGrid(const GridType& knots)
        {
          bKnots=knots;
        }

        ///empty reset
        void reset()
        {
        }

        inline void setGrid(real_type xi, real_type xf, 
            real_type yi, real_type yf, real_type zi, real_type zf, 
            int nx, int ny, int nz, 
            bool interp=true, bool periodic=true,bool openend=true)
        {
          bKnots.setGrid(xi,xf,yi,yf,zi,zf,nx,ny,nz,interp,periodic,openend);
        }

        /** add a orbital
         * @param i index of the orbital
         * @param data input data
         * @param curP interpolated data
         */
        void add(int i, const StorageType& data, StorageType* curP)
        {
          if(i<OffSet) OffSet=i;
          DeleteP.push_back(false);
          P.push_back(curP);
          bKnots.Init(data,*curP);
        }

        /** add a orbital
         * @param i index of the orbital
         * @param data input data
         *
         * New interpolated data is created and will be deleted by the constructor.
         */
        void add(int i, const StorageType& data)
        {
          if(i<OffSet) OffSet=i;
          StorageType *curP=new StorageType;
          DeleteP.push_back(true);
          P.push_back(curP);
          bKnots.Init(data,*curP);
        }

        template<typename PV>
        inline void evaluate(const TinyVector<real_type,3>& r, PV& vals) 
        {
          bKnots.Find(r[0],r[1],r[2]);
          for(int m=0, j=OffSet;m<P.size(); m++,j++)
          {
            vals[j]=bKnots.evaluate(*P[m]);
          }
        }

        template<typename PV, typename GV>
        inline void
          evaluate(const TinyVector<real_type,3>& r, PV& vals, GV& grads, PV& laps)
          {
            bKnots.FindAll(r[0],r[1],r[2]);
            for(int m=0,j=OffSet;m<P.size(); m++,j++)
            {
              vals[j]=bKnots.evaluate(*P[m],grads[j],laps[j]);
            }
          }

        template<typename PM, typename GM>
        inline void
          evaluate(const TinyVector<real_type,3>& r, int i, PM& vals, GM& grads, PM& laps)
          {
            bKnots.FindAll(r[0],r[1],r[2]);
            for(int m=0,j=OffSet;m<P.size(); m++,j++)
            {
              vals(j,i)=bKnots.evaluate(*P[m],grads(i,j),laps(i,j));
            }
          }

      private:
        int OffSet;
        GridType bKnots;
        std::vector<bool> DeleteP;
        std::vector<StorageType*> P;
    };
}

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/
#endif

