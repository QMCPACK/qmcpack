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
    
    
#ifndef QMCPLUSPLUS_TRICUBIC_B_SPLINE_TWISTSET_H
#define QMCPLUSPLUS_TRICUBIC_B_SPLINE_TWISTSET_H

#include "Numerics/TricubicBsplineGrid.h"
#include <map>

namespace qmcplusplus
{

/** A group of bspline functions stored in a std::map<int,StorageType*> for a twist angle
 */
template<typename T>
class TricubicBsplineTwistSet: public TricubicBsplineTraits<T>
{
public:
  typedef typename TricubicBsplineTraits<T>::real_type   real_type;
  typedef typename TricubicBsplineTraits<T>::value_type  value_type;
  typedef typename TricubicBsplineTraits<T>::PosType     PosType;
  typedef typename TricubicBsplineTraits<T>::GridType    GridType;
  typedef typename TricubicBsplineTraits<T>::StorageType StorageType;
  typedef typename std::map<int,const StorageType*>::iterator  IteratorType;
  typedef TinyVector<value_type,3> GradType;

  //going to use GGt=dot(Lattice.G,transpose(Lattice.G))
  using TricubicBsplineTraits<T>::GGt;
  //going to use Lattice
  using TricubicBsplineTraits<T>::Lattice;

  /** default constructure
   */
  TricubicBsplineTwistSet()
  {
  }

  ~TricubicBsplineTwistSet()
  {
  }

  /** Set the Cartesian Twist angle
   * @param tangle Twist angle in the Cartesian Corrdinate
   */
  inline void setTwistAngle(const PosType& tangle)
  {
    TwistAngle=tangle;
    mK2=-dot(tangle,tangle);
  }

  inline void setGrid(const GridType& knots)
  {
    bKnots=knots;
  }

  ///empty reset
  void resetParameters(VarRegistry<real_type>& vlist)
  {
    ///DO NOTHING FOR NOW
  }

  inline void setGrid(
    real_type xi, real_type xf,
    real_type yi, real_type yf,
    real_type zi, real_type zf,
    int nx, int ny, int nz,
    bool interp=true, bool periodic=true,bool openend=true)
  {
    //bKnots.setGrid(xi,xf,yi,yf,zi,zf,nx,ny,nz,interp,periodic,openend);
    bKnots.setGrid(0.0,1.0,0.0,1.0,0.0,1.0,nx,ny,nz,interp,periodic,openend);
  }

  /** add a orbital
   * @param i index of the orbital
   * @param data input data
   * @param curP interpolated data
   */
  void add(int i, const PosType& c, const StorageType& data, StorageType* curP)
  {
    IteratorType pit(P.find(i));
    if(pit == P.end())
    {
      bKnots.Init(data,*curP);
      P[i]=curP;
    }
  }

  void add(int i, const PosType& c, StorageType* curP)
  {
    IteratorType pit(P.find(i));
    if(pit == P.end())
    {
      P[i]=curP;
    }
  }

  template<typename PV>
  inline void evaluate(const PosType& r, PV& vals)
  {
    PosType ru(Lattice.toUnit(r));
    bKnots.Find(ru[0],ru[1],ru[2]);
    real_type phi(dot(TwistAngle,r));
    value_type phase(std::cos(phi),std::sin(phi));
    IteratorType pit(P.begin()), pit_end(P.end());
    while(pit != pit_end)
    {
      vals[(*pit).first]=phase*bKnots.evaluate(*((*pit).second));
      ++pit;
    }
  }

  template<typename PV, typename GV>
  inline void
  evaluate(const PosType& r, PV& vals, GV& grads, PV& laps)
  {
    PosType ru(Lattice.toUnit(r));
    bKnots.FindAll(ru[0],ru[1],ru[2]);
    real_type phi(dot(TwistAngle,r));
    real_type c=std::cos(phi),s=std::sin(phi);
    value_type phase(c,s);
    GradType dk(value_type(-TwistAngle[0]*s,TwistAngle[0]*c),
                value_type(-TwistAngle[1]*s,TwistAngle[1]*c),
                value_type(-TwistAngle[2]*s,TwistAngle[2]*c));
    TinyVector<value_type,3> gu;
    Tensor<value_type,3> hess;
    IteratorType pit(P.begin()), pit_end(P.end());
    while(pit != pit_end)
    {
      int j((*pit).first);
      value_type v=bKnots.evaluate(*((*pit).second),gu,hess);
      GradType g=dot(Lattice.G,gu);
      value_type l=trace(hess,GGt);
      vals[j]=phase*v;
      grads[j]=v*dk+phase*g;
      laps[j]=phase*(mK2*v+l)+2.0*dot(dk,g);
      ++pit;
    }
  }

  template<typename PM, typename GM>
  inline void
  evaluate(const PosType& r, int i, PM& vals, GM& grads, PM& laps)
  {
    PosType ru(Lattice.toUnit(r));
    bKnots.FindAll(ru[0],ru[1],ru[2]);
    real_type phi(dot(TwistAngle,r));
    real_type c=std::cos(phi),s=std::sin(phi);
    value_type phase(c,s);
    GradType dk(value_type(-TwistAngle[0]*s,TwistAngle[0]*c),
                value_type(-TwistAngle[1]*s,TwistAngle[1]*c),
                value_type(-TwistAngle[2]*s,TwistAngle[2]*c));
    TinyVector<value_type,3> gu;
    Tensor<value_type,3> hess;
    IteratorType pit(P.begin()), pit_end(P.end());
    while(pit != pit_end)
    {
      int j((*pit).first);
      value_type v=bKnots.evaluate(*((*pit).second),gu,hess);
      GradType g=dot(Lattice.G,gu);
      value_type l=trace(hess,GGt);
      vals(j,i)=phase*v;
      grads(i,j)=v*dk+phase*g;
      laps(i,j)=phase*(mK2*v+l)+2.0*dot(dk,g);
      ++pit;
    }
  }

private:
  real_type mK2;
  PosType TwistAngle;
  GridType bKnots;
  std::map<int,const StorageType*> P;
};

//  /** A group of bspline functions stored in a std::map<int,StorageType*>
//   */
//  template<typename T>
//    class TricubicBsplineGCSet
//    {
//      public:
//        typedef typename TricubicBsplineTraits<T>::real_type   real_type;
//        typedef typename TricubicBsplineTraits<T>::value_type  value_type;
//        typedef typename TricubicBsplineTraits<T>::PosType     PosType;
//        typedef typename TricubicBsplineTraits<T>::GridType    GridType;
//        typedef typename TricubicBsplineTraits<T>::StorageType StorageType;
//        typedef typename std::map<int,const StorageType*>::iterator  IteratorType;
//        typedef TinyVector<value_type,3> GradType;
//        typedef Tensorr<value_type,3>    TensorType;
//
//        /** default constructure
//         */
//        TricubicBsplineGCSet()
//        {
//        }
//
//        ~TricubicBsplineGCSet() {
//        }
//
//        inline void setGrid(const GridType& knots)
//        {
//          bKnots=knots;
//        }
//
//        ///empty reset
//        void resetParameters(VarRegistry<real_type>& vlist)
//        {
//          ///DO NOTHING FOR NOW
//        }
//
//        /** set tricubic grid
//         *
//         * The grid is set to [0,1]x[0,1]x[0,2]
//         */
//        template<typename LT>
//        inline void setGrid(const LT& lat, int nx, int ny, int nz,
//            bool interp=true, bool periodic=true,bool openend=true)
//        {
//          Lattice=lat;
//          J = Lattice.G;//transpose?
//          bKnots.setGrid(xi,1.0,yi,1.0,zi,1.0,nx,ny,nz,interp,periodic,openend);
//        }
//
//        /** add a orbital
//         * @param i index of the orbital
//         * @param data input data
//         * @param curP interpolated data
//         */
//        void add(int i, const StorageType& data, StorageType* curP)
//        {
//          IteratorType pit(P.find(i));
//          if(pit == P.end())
//          {
//            bKnots.Init(data,*curP);
//            P[i]=curP;
//          }
//        }
//
//        void add(int i, StorageType* curP)
//        {
//          IteratorType pit(P.find(i));
//          if(pit == P.end())
//          {
//            P[i]=curP;
//          }
//        }
//
//        template<typename PV>
//        inline void evaluate(const PosType& r, PV& vals)
//        {
//          PosType ru=Lattice.toUnit(r);
//          bKnots.Find(ru[0],ru[1],ru[2]);
//          IteratorType pit(P.begin()), pit_end(P.end());
//          while(pit != pit_end)
//          {
//            vals[(*pit).first]=bKnots.evaluate(*((*pit).second));
//            ++pit;
//          }
//        }
//
//        template<typename PV, typename GV>
//        inline void
//          evaluate(const PosType& r, PV& vals, GV& grads, PV& laps)
//          {
//            PosType ru=Lattice.toUnit(r);
//            bKnots.FindAll(ru[0],ru[1],ru[2]);
//
//            IteratorType pit(P.begin()), pit_end(P.end());
//            while(pit != pit_end)
//            {
//              int j((*pit).first);
//              GradType gr;
//              TensorType t3;
//              vals[j]=bKnots.evaluate(*((*pit).second),gr,t3);
//              grads[j]=dot(J,gr);//POTENTIAL BUG WITH COMPLEX
//              laps[j]=J(0,0)*(t3(0,0)*J(0,0)+t3(0,1)*J(0,1)+t3(0,2)*J(0,2))
//                     +J(0,1)*(t3(1,0)*J(0,0)+t3(1,1)*J(0,1)+t3(1,2)*J(0,2))
//                     +J(0,2)*(t3(2,0)*J(0,0)+t3(2,1)*J(0,1)+t3(2,2)*J(0,2))
//                     +J(1,0)*(t3(0,0)*J(0,0)+t3(0,1)*J(0,1)+t3(0,2)*J(0,2))
//                     +J(1,1)*(t3(1,0)*J(0,0)+t3(1,1)*J(0,1)+t3(1,2)*J(0,2))
//                     +J(1,2)*(t3(2,0)*J(0,0)+t3(2,1)*J(0,1)+t3(2,2)*J(0,2))
//                     +J(2,0)*(t3(0,0)*J(0,0)+t3(0,1)*J(0,1)+t3(0,2)*J(0,2))
//                     +J(2,1)*(t3(1,0)*J(0,0)+t3(1,1)*J(0,1)+t3(1,2)*J(0,2))
//                     +J(2,2)*(t3(2,0)*J(0,0)+t3(2,1)*J(0,1)+t3(2,2)*J(0,2));
//              //vals[j]=bKnots.evaluate(*((*pit).second),grads[j],laps[j]);
//              ++pit;
//            }
//          }
//
//        template<typename PM, typename GM>
//        inline void
//          evaluate(const PosType& r, int i, PM& vals, GM& grads, PM& laps)
//          {
//            PosType ru=Lattice.toUnit(r);
//            bKnots.FindAll(ru[0],ru[1],ru[2]);
//            IteratorType pit(P.begin()), pit_end(P.end());
//            while(pit != pit_end)
//            {
//              int j((*pit).first);
//              GradType gr;
//              TensorType t3;
//              vals(j,i)=bKnots.evaluate(*((*pit).second),gr,t3);
//              grads(i,j)=dot(J,gr);
//              laps(i,j)
//                    = J(0,0)*(t3(0,0)*J(0,0)+t3(0,1)*J(0,1)+t3(0,2)*J(0,2))
//                     +J(0,1)*(t3(1,0)*J(0,0)+t3(1,1)*J(0,1)+t3(1,2)*J(0,2))
//                     +J(0,2)*(t3(2,0)*J(0,0)+t3(2,1)*J(0,1)+t3(2,2)*J(0,2))
//                     +J(1,0)*(t3(0,0)*J(0,0)+t3(0,1)*J(0,1)+t3(0,2)*J(0,2))
//                     +J(1,1)*(t3(1,0)*J(0,0)+t3(1,1)*J(0,1)+t3(1,2)*J(0,2))
//                     +J(1,2)*(t3(2,0)*J(0,0)+t3(2,1)*J(0,1)+t3(2,2)*J(0,2))
//                     +J(2,0)*(t3(0,0)*J(0,0)+t3(0,1)*J(0,1)+t3(0,2)*J(0,2))
//                     +J(2,1)*(t3(1,0)*J(0,0)+t3(1,1)*J(0,1)+t3(1,2)*J(0,2))
//                     +J(2,2)*(t3(2,0)*J(0,0)+t3(2,1)*J(0,1)+t3(2,2)*J(0,2));
//              //vals(j,i)=bKnots.evaluate(*((*pit).second),grads(i,j),laps(i,j));
//              ++pit;
//            }
//          }
//
//      private:
//        GridType bKnots;
//        PosType J;
//        LatticeType Lattice;
//        std::map<int,const StorageType*> P;
//    };
}
#endif
