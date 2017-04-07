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
    
    



#ifndef QMCPLUSPLUS_TRICUBICSPLINETEMPLATE_H
#define QMCPLUSPLUS_TRICUBICSPLINETEMPLATE_H
#include "Numerics/XYZCubicGrid.h"
#include "Numerics/OneDimCubicSpline.h"

namespace qmcplusplus
{
/** Tri-cubic Splines with periodic boundary conditions and fixed first derivatives.
 *
 * Adapting TriCubicSpline implemented by K. Esler and D. Das.
 * Use stl containers
 */
template<class T, class Tg=T>
class TriCubicSplineT
{

public:

  typedef T                             value_type;
  typedef Tg                            point_type;
  typedef XYZCubicGrid<T,Tg>            GridType;
  typedef TriCubicSplineT<T,Tg>         ThisType;
  typedef typename GridType::KnotType   KnotType;
  typedef typename GridType::Grid1DType Grid1DType;

  //true, if this function evaluates grid-related properties
  bool GridManager;
  //Hide these later
  bool UpToDate;
  int nX, nY, nZ;
  int n001,n010,n011,n100,n101,n110,n111;
  GridType* m_grid;
  std::vector<KnotType> F;

  /// constructor
  TriCubicSplineT(GridType* agrid):
    GridManager(true),UpToDate(false), nX(0), nY(0), nZ(0), m_grid(agrid)
  {
    if(m_grid)
    {
      nX = m_grid->nX;
      nY = m_grid->nY;
      nZ = m_grid->nZ;
      F.resize(nX*nY*nZ);
      n100=nY*nZ;
      n001=1;
      n010=nZ;
      n011=nZ+1;
      n101=n100+n001;
      n110=n100+n010;
      n111=n100+n011;
    }
  }

  inline void setGridManager(bool willmanage)
  {
    GridManager=willmanage;
  }

  inline void reset(bool periodic=true)
  {
    if(m_grid)
    {
      m_grid->setBC(periodic);
      if(!UpToDate)
      {
        UpdateX(0, 1, periodic);  // Do dF/dx
        UpdateY(0, 2, periodic);  // Do dF/dy
        UpdateZ(0, 3, periodic);  // Do dF/dz
        UpdateY(1, 4, periodic);  // Do d2F/dxdy
        UpdateZ(1, 5, periodic);  // Do d2F/dxdz
        UpdateZ(2, 6, periodic);  // Do d2F/dydz
        UpdateZ(4, 7, periodic);  // Do d3F/dxdydz
        UpToDate=true;
      }
    }
  }

  inline T* data()
  {
    return &(F[0][0]);
  }
  inline const T* data() const
  {
    return &(F[0][0]);
  }

  inline T operator()(int n) const
  {
    return F[n][0];
  }

  inline T& operator()(int n)
  {
    return F[n][0];
  }

  inline T operator()(int i, int j, int k) const
  {
    return F[index(i,j,k)][0];
  }

  inline T& operator()(int i, int j, int k)
  {
    return F[index(i,j,k)][0];
  }

  inline int index(int i, int j, int k)
  {
    return k+nZ*(j+nY*i);
  }

  template<class IT>
  inline void reset(IT first, IT last, bool periodic=true)
  {
    typename std::vector<KnotType>::iterator it(F.begin());
    while(first != last)
    {
      (*it)[0]=*first++;
      ++it;
    }
    reset(periodic);
  }

  //template<class PV>
  //inline void setgrid(const PV& r) {
  //  m_grid->locate(r[0],r[1],r[2]);
  //}

  inline T evaluate(const TinyVector<Tg,3>& r)
  {
    if(GridManager)
    {
      m_grid->locate(r[0],r[1],r[2],false);
      //m_grid->update(false);
    }
    if(m_grid->Loc<0)
    {
      return 1e-20;
    }
    //m_grid->update();
    int cur(m_grid->Loc);
    return
      m_grid->evaluate(F[cur]    , F[cur+n001], F[cur+n010], F[cur+n011],
                       F[cur+n100], F[cur+n101], F[cur+n110], F[cur+n111]);
  }

  inline T evaluate(const TinyVector<Tg,3>& r, TinyVector<T,3>& gradf, T& lapf)
  {
    if(GridManager)
    {
      m_grid->locate(r[0],r[1],r[2],true);
      //m_grid->update(true);
    }
    if(m_grid->Loc<0)
    {
      gradf[0] = 1e-40;
      gradf[1] = 1e-40;
      gradf[2] = 1e-40;
      lapf = 1e-40;
      return 1e-20;
    }
    //m_grid->updateAll();
    int cur(m_grid->Loc);
    m_grid->evaluateAll(F[cur]    ,  F[cur+n001], F[cur+n010], F[cur+n011],
                        F[cur+n100], F[cur+n101], F[cur+n110], F[cur+n111]);
    gradf[0]=m_grid->gradfX;
    gradf[1]=m_grid->gradfY;
    gradf[2]=m_grid->gradfZ;
    lapf=m_grid->lapf;
    return m_grid->val;
  }

  void resetParameters(VarRegistry<point_type>& vlist)
  {
    ///DO NOTHING FOR NOW
  }

private:


  // dim:     Dimension to calculate derivative w.r.t
  // source:  Function to differentiate
  // dest:    where to put result
  void UpdateX (int source, int target, bool periodic);
  void UpdateY (int source, int target, bool periodic);
  void UpdateZ (int source, int target, bool periodic);
};

template<class T, class Tg>
void TriCubicSplineT<T,Tg>::UpdateX(int source, int target, bool periodic)
{
  OneDimGridFunctor<T,Tg>* temp=0;
  if(periodic)
  {
    temp = new OneDimCubicSplinePBC<T,Tg>(m_grid->gridX);
  }
  else
  {
    temp = new OneDimCubicSplineFirst<T,Tg>(m_grid->gridX);
  }
  int nyz=nY*nZ;
  // Loop over all y and z
  for(int iy = 0; iy < nY; iy++)
  {
    for(int iz = 0; iz < nZ; iz++)
    {
      int first(index(0,iy,iz));
      T* restrict s = temp->data();
      for(int ix=0,cur=first; ix< nX; ix++,cur+=nyz,s++)
        *s = F[cur][source];
      temp->spline();
      const T* restrict sfit = temp->data(1);
      for(int ix=0,cur=first; ix< nX; ix++,cur+=nyz)
        F[cur][target] = *sfit++;
    }
  }
  delete temp;
}

template<class T, class Tg>
void TriCubicSplineT<T,Tg>::UpdateY(int source, int target, bool periodic)
{
  OneDimGridFunctor<T,Tg>* temp=0;
  if(periodic)
  {
    temp = new OneDimCubicSplinePBC<T,Tg>(m_grid->gridY);
  }
  else
  {
    temp = new OneDimCubicSplineFirst<T,Tg>(m_grid->gridY);
  }
  // Loop over all x and z
  for(int ix = 0; ix < nX; ix++)
  {
    for(int iz = 0; iz < nZ; iz++)
    {
      int first(index(ix,0,iz));
      T* restrict s = temp->data();
      for(int iy=0,cur=first; iy< nY; iy++,cur+=nZ,s++)
        *s = F[cur][source];
      temp->spline();
      const T* restrict sfit = temp->data(1);
      for(int iy=0,cur=first; iy< nY; iy++,cur+=nZ)
        F[cur][target] = *sfit++;
    }
  }
  delete temp;
}

template<class T, class Tg>
void TriCubicSplineT<T,Tg>::UpdateZ(int source, int target, bool periodic)
{
  OneDimGridFunctor<T,Tg>* temp=0;
  if(periodic)
  {
    temp = new OneDimCubicSplinePBC<T,Tg>(m_grid->gridZ);
  }
  else
  {
    temp = new OneDimCubicSplineFirst<T,Tg>(m_grid->gridZ);
  }
  // Loop over all x and z
  for(int ix = 0; ix < nX; ix++)
  {
    for(int iy = 0; iy < nY; iy++)
    {
      int first(index(ix,iy,0));
      T* restrict s = temp->data();
      for(int iz=0,cur=first; iz< nZ; iz++,cur++,s++)
        *s = F[cur][source];
      temp->spline();
      const T* restrict sfit = temp->data(1);
      for(int iz=0,cur=first; iz< nZ; iz++,cur++)
        F[cur][target] = *sfit++;
    }
  }
  delete temp;
}

}
#endif
