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

#include "Numerics/TricubicBsplineGrid.h"

namespace qmcplusplus {

  template<typename T>
    class TricubicBspline
    {
      public:
        typedef OrbitalTraits<T>::real_type real_type;
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
            return bKnots.evaluate(P,grad,lapf);
          }

      private:
        //Grid
        GridType bKnots;
        // The control points
        Array<T,3> P;
    };

  template<typename T>
    class TricubicBsplineSet
    {
      public:
        typedef OrbitalTraits<T>::real_type real_type;
        typedef TricubicBsplineGrid<T> GridType;

        TricubicBsplineSet() { }
        ~TricubicBsplineSet() { 
          delete_iter(P.begin(), P.end());
        }

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

        void add(int i, const Array<T,3>& data)
        {
          int offset=P.size();
          Array<T,3> *curP=new Array<T,3>;
          P.push_back(curP);
          Pid.push_back(i);
          bKnots.Init(data,*curP);
        }

        template<typename PV>
        inline void evaluate(const TinyVector<real_type,3>& r, PV& vals) 
        {
          bKnots.Find(r[0],r[1],r[2]);
          for(int i=0;i<Pid.size(); i++)
          {
            vals[Pid[i]]=bKnots.evaluate(*P[i]);
          }
        }

        template<typename PV, typename GV>
        inline void
          evaluate(const TinyVector<real_type,3>& r, PV& vals, GV& grads, PV& laps)
          {
            bKnots.FindAll(r[0],r[1],r[2]);
            for(int i=0;i<Pid.size(); i++)
            {
              int j=Pid[i];
              vals[j]=bKnots.evaluate(*P[i],grads[j],laps[j]);
            }
          }

      private:
        GridType bKnots;
        std::vector<Array<T,3>*> P;
    };
}

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/
#endif

