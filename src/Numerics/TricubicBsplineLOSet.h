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
    
    



#ifndef TRICUBIC_B_SPLINE_LOCALIZEDORBITAL_SET_H
#define TRICUBIC_B_SPLINE_LOCALIZEDORBITAL_SET_H

#include "Numerics/TricubicBsplineGrid.h"
#include <map>

namespace qmcplusplus
{

/** A group of bspline functions stored in a std::map<int,StorageType*>
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
  std::vector<PosType> Centers;

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
