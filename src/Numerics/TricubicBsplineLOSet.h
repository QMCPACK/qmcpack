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
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef TRICUBIC_B_SPLINE_LOCALIZEDORBITAL_SET_H
#define TRICUBIC_B_SPLINE_LOCALIZEDORBITAL_SET_H

#include "Numerics/TricubicBsplineGrid.h"
#include <map>

namespace qmcplusplus
{

/** A group of bspline functions stored in a map<int,StorageType*>
 */
template<class T>
class TricubicBsplineLOSet: public TricubicBsplineTraits<T>
{
public:
  typedef typename TricubicBsplineTraits<T>::real_type   real_type;
  typedef typename TricubicBsplineTraits<T>::value_type  value_type;
  typedef typename TricubicBsplineTraits<T>::PosType     PosType;
  typedef typename TricubicBsplineTraits<T>::GridType    GridType;
  typedef typename TricubicBsplineTraits<T>::StorageType StorageType;
  typedef typename std::map<int,const StorageType*>::iterator  IteratorType;

  using TricubicBsplineTraits<T>::Rcut2;
  vector<PosType> Centers;

  /** default constructure
   *
   * Set Rcut2 to a large number so that everything counts
   */
  TricubicBsplineLOSet()
  {
    Rcut2=1e6;
  }

  ~TricubicBsplineLOSet()
  {
  }

  inline void setTwistAngle(const PosType& tangle)
  {
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

  inline void setGrid(real_type xi, real_type xf,
                      real_type yi, real_type yf, real_type zi, real_type zf,
                      int nx, int ny, int nz,
                      bool interp=true, bool periodic=true,bool openend=true)
  {
    bKnots.setGrid(xi,xf,yi,yf,zi,zf,nx,ny,nz,interp,periodic,openend);
    //ParticleSet::Tensor_t sc(xf-xi,0.0,0.0,0.0,yf-yi,0.0,0.0,0.0,zf-zi);
    //Lattice.set(sc);
  }

  /** add a orbital
   * @param i index of the orbital
   * @param data input data
   * @param curP interpolated data
   */
  void add(int i, const PosType& c, const StorageType& data, StorageType* curP)
  {
    bKnots.Init(data,*curP);
    Centers.push_back(c);
    P.push_back(curP);
  }

  void add(int i, const PosType& c, StorageType* curP)
  {
    //Already exists
    if(i<Centers.size())
      return;
    Centers.push_back(c);
    P.push_back(curP);
  }

  template<typename PV>
  inline void evaluate(const PosType& r, PV& vals)
  {
    bKnots.Find(r[0],r[1],r[2]);
    for(int j=0; j<Centers.size(); j++)
    {
      if(bKnots.getSep2(r[0]-Centers[j][0],r[1]-Centers[j][1],r[2]-Centers[j][2])>Rcut2)
        vals[j]=0.0;//numeric_limits<T>::epsilon();
      else
        vals[j]=bKnots.evaluate(*P[j]);
    }
  }

  template<typename PV, typename GV>
  inline void
  evaluate(const PosType& r, PV& vals, GV& grads, PV& laps)
  {
    bKnots.FindAll(r[0],r[1],r[2]);
    for(int j=0; j<Centers.size(); j++)
    {
      if(bKnots.getSep2(r[0]-Centers[j][0],r[1]-Centers[j][1],r[2]-Centers[j][2])>Rcut2)
      {
        vals[j]=0.0;//numeric_limits<T>::epsilon();
        grads[j]=0.0;
        laps[j]=0.0;
      }
      else
        vals[j]=bKnots.evaluate(*P[j],grads[j],laps[j]);
    }
  }

  template<typename PM, typename GM>
  inline void
  evaluate(const PosType& r, int i, PM& vals, GM& grads, PM& laps)
  {
    bKnots.FindAll(r[0],r[1],r[2]);
    for(int j=0; j<Centers.size(); j++)
    {
      if(bKnots.getSep2(r[0]-Centers[j][0],r[1]-Centers[j][1],r[2]-Centers[j][2])>Rcut2)
      {
        vals(j,i)=0.0; //numeric_limits<T>::epsilon();
        grads(i,j)=0.0;
        laps(i,j)=0.0;
      }
      else
        vals(j,i)=bKnots.evaluate(*P[j],grads(i,j),laps(i,j));
    }
  }

private:
  GridType bKnots;
  std::vector<const StorageType*> P;
};
}
#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1863 $   $Date: 2007-04-01 12:36:53 -0500 (Sun, 01 Apr 2007) $
 * $Id: TricubicBsplineLOSet.h 1863 2007-04-01 17:36:53Z jnkim $
 ***************************************************************************/
