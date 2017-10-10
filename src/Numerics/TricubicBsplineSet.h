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
    
    



#ifndef TRICUBIC_B_SPLINE_SET_H
#define TRICUBIC_B_SPLINE_SET_H

#include "Numerics/TricubicBsplineGrid.h"
#include <map>

namespace qmcplusplus
{

template<typename T>
class TricubicBspline: public TricubicBsplineTraits<T>
{
public:
  typedef typename TricubicBsplineTraits<T>::real_type   real_type;
  typedef typename TricubicBsplineTraits<T>::value_type  value_type;
  typedef typename TricubicBsplineTraits<T>::PosType     PosType;
  typedef typename TricubicBsplineTraits<T>::GridType    GridType;
  typedef typename TricubicBsplineTraits<T>::StorageType StorageType;

  TricubicBspline() {}

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


/** A group of bspline functions stored in a std::vector<StorageType*>
 *
 * Assume a linear order of the bspline sets.
 */
template<typename T>
class TricubicBsplineVector: public TricubicBsplineTraits<T>
{
public:
  typedef typename TricubicBsplineTraits<T>::real_type   real_type;
  typedef typename TricubicBsplineTraits<T>::value_type  value_type;
  typedef typename TricubicBsplineTraits<T>::PosType     PosType;
  typedef typename TricubicBsplineTraits<T>::GridType    GridType;
  typedef typename TricubicBsplineTraits<T>::StorageType StorageType;

  /** default constructure
   *
   * For memory efficiency, reserve DeleteP and P
   * OffSet is set to 1000000. Safe until we can do 1000000 orbitals.
   */
  TricubicBsplineVector():OffSet(1000000)
  {
    DeleteP.reserve(1024);
    P.reserve(1024);
  }

  ~TricubicBsplineVector()
  {
    for(int i=0; i<DeleteP.size(); i++)
    {
      if(DeleteP[i])
        delete P[i];
    }
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
  }

  /** add a orbital
   * @param i index of the orbital
   * @param c center dummy argument to work with truncated case
   * @param data input data
   * @param curP interpolated data
   */
  void add(int i, const PosType& c, const StorageType& data, StorageType* curP)
  {
    if(i<OffSet)
      OffSet=i;
    DeleteP.push_back(false);
    P.push_back(curP);
    bKnots.Init(data,*curP);
  }

  /** add a orbital
   * @param i index of the orbital
   * @param c center dummy argument to work with truncated case
   * @param data input data
   *
   * New interpolated data is created and will be deleted by the constructor.
   */
  void add(int i,const PosType& c, const StorageType& data)
  {
    if(i<OffSet)
      OffSet=i;
    StorageType *curP=new StorageType;
    DeleteP.push_back(true);
    P.push_back(curP);
    bKnots.Init(data,*curP);
  }

  inline value_type
  evaluate(int iorb, const PosType& r)
  {
    bKnots.Find(r[0],r[1],r[2]);
    return bKnots.evaluate(*P[iorb]);
  }

  template<typename PV>
  inline void evaluate(const PosType& r, PV& vals)
  {
    bKnots.Find(r[0],r[1],r[2]);
    for(int m=0, j=OffSet; m<P.size(); m++,j++)
    {
      vals[j]=bKnots.evaluate(*P[m]);
    }
  }

  template<typename PV, typename GV>
  inline void
  evaluate(const PosType& r, PV& vals, GV& grads, PV& laps)
  {
    bKnots.FindAll(r[0],r[1],r[2]);
    for(int m=0,j=OffSet; m<P.size(); m++,j++)
    {
      vals[j]=bKnots.evaluate(*P[m],grads[j],laps[j]);
    }
  }

  template<typename PM, typename GM>
  inline void
  evaluate(const PosType& r, int i, PM& vals, GM& grads, PM& laps)
  {
    bKnots.FindAll(r[0],r[1],r[2]);
    for(int m=0,j=OffSet; m<P.size(); m++,j++)
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

/** A group of bspline functions stored in a std::map<int,StorageType*>
 */
template<typename T>
class TricubicBsplineSet: public TricubicBsplineTraits<T>
{
public:
  typedef typename TricubicBsplineTraits<T>::real_type   real_type;
  typedef typename TricubicBsplineTraits<T>::value_type  value_type;
  typedef typename TricubicBsplineTraits<T>::PosType     PosType;
  typedef typename TricubicBsplineTraits<T>::GridType    GridType;
  typedef typename TricubicBsplineTraits<T>::StorageType StorageType;
  typedef typename std::map<int,const StorageType*>::iterator  IteratorType;

  /** default constructure
   */
  TricubicBsplineSet()
  {
  }

  ~TricubicBsplineSet()
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

  void add(int i,const PosType& c,  StorageType* curP)
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
    bKnots.Find(r[0],r[1],r[2]);
    IteratorType pit(P.begin()), pit_end(P.end());
    while(pit != pit_end)
    {
      vals[(*pit).first]=bKnots.evaluate(*((*pit).second));
      ++pit;
    }
  }

  template<typename PV, typename GV>
  inline void
  evaluate(const PosType& r, PV& vals, GV& grads, PV& laps)
  {
    bKnots.FindAll(r[0],r[1],r[2]);
    IteratorType pit(P.begin()), pit_end(P.end());
    while(pit != pit_end)
    {
      int j((*pit).first);
      vals[j]=bKnots.evaluate(*((*pit).second),grads[j],laps[j]);
      ++pit;
    }
  }

  template<typename PM, typename GM>
  inline void
  evaluate(const PosType& r, int i, PM& vals, GM& grads, PM& laps)
  {
    bKnots.FindAll(r[0],r[1],r[2]);
    IteratorType pit(P.begin()), pit_end(P.end());
    while(pit != pit_end)
    {
      int j((*pit).first);
      vals(j,i)=bKnots.evaluate(*((*pit).second),grads(i,j),laps(i,j));
      ++pit;
    }
  }

private:
  GridType bKnots;
  std::map<int,const StorageType*> P;
};
}
#endif
