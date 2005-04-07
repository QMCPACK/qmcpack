#ifndef OHMMS_QMC_TRICUBICSPLINETEMPLATE_H
#define OHMMS_QMC_TRICUBICSPLINETEMPLATE_H
#include "Numerics/XYZCubicGrid.h"
#include "Numerics/OneDimCubicSpline.h"
/** Tri-cubic Splines with periodic boundary conditions and fixed first derivatives.
 *
 * Adapting TriCubicSpline implemented by K. Esler and D. Das.
 * Use stl containers
 */
template<class T>
class TriCubicSplineT {

public: 

  typedef XYZCubicGrid<T>               GridType;
  typedef typename GridType::KnotType   KnotType;
  typedef typename GridType::Grid1DType Grid1DType;
  typedef T value_type;

  /// constructor
  TriCubicSplineT(GridType* agrid): 
    UpToDate(false), nX(0), nY(0), nZ(0), m_grid(agrid)
  {
    if(m_grid) {
      nX = m_grid->nX;
      nY = m_grid->nY;
      nZ = m_grid->nZ;
      F.resize(nX*nY*nZ);
      n100=nY*nZ; 
      n001=1; n010=nZ;  n011=nZ+1;
      n101=n100+n001; n110=n100+n010; n111=n100+n011;
    }
  }

  inline void reset() { 
    if(m_grid) {
      if(!UpToDate) {
        UpdateX(0, 1);  // Do dF/dx
        UpdateY(0, 2);  // Do dF/dy
        UpdateZ(0, 3);  // Do dF/dz
        UpdateY(1, 4);  // Do d2F/dxdy
        UpdateZ(1, 5);  // Do d2F/dxdz
        UpdateZ(2, 6);  // Do d2F/dydz
        UpdateZ(4, 7);  // Do d3F/dxdydz
        UpToDate=true;
      }
    } 
  }

  inline T operator()(int i, int j, int k) const
  {
    return F[index(i,j,k)][0];
  }
  
  inline T& operator()(int i, int j, int k) 
  {
    return F[index(i,j,k)][0];
  }

  inline int index(int i, int j, int k) {
    return k+nZ*(j+nY*i);
  }

  template<class PV>
  inline void setgrid(const PV& r) {
    m_grid->locate(r[0],r[1],r[2]);
  }

  template<class PV>
  inline T evaluate(const PV& r){
    if(m_grid->Loc<0) {
      return 1e-20;
    }
    m_grid->update();
    int cur(m_grid->Loc);
    m_grid->evaluate(F[cur]    , F[cur+n001], F[cur+n010], F[cur+n011],
                    F[cur+n100], F[cur+n101], F[cur+n110], F[cur+n111]);
    return m_grid->val;
  }

  template<class PV>
  inline T evaluate(const PV& r, PV& gradf, T& lapf){
    if(m_grid->Loc<0){
      gradf[0] = 1e-40; gradf[1] = 1e-40; gradf[2] = 1e-40;
      lapf = 1e-40;
      return 1e-20;
    }
    m_grid->updateAll();
    int cur(m_grid->Loc);
    m_grid->evaluateAll(F[cur]    , F[cur+n001], F[cur+n010], F[cur+n011],
                        F[cur+n100], F[cur+n101], F[cur+n110], F[cur+n111]);
        
    gradf[0]=m_grid->gradfX; gradf[1]=m_grid->gradfY; gradf[2]=m_grid->gradfZ;
    lapf=m_grid->lapf;
    return m_grid->val;
  }

private:

  bool UpToDate; 
  int nX, nY, nZ;
  int n001,n010,n011,n100,n101,n110,n111;

  GridType* m_grid;
  std::vector<KnotType> F;

  // dim:     Dimension to calculate derivative w.r.t
  // source:  Function to differentiate
  // dest:    where to put result
  void UpdateX (int source, int target, bool periodic=true);
  void UpdateY (int source, int target, bool periodic=true);
  void UpdateZ (int source, int target, bool periodic=true);
};

template<class T>
void TriCubicSplineT<T>::UpdateX(int source, int target, bool periodic){
  
  typename OneDimGridFunctor<T>* temp=0;
  if(periodic) {
    temp = new OneDimCubicSplinePBC<T>(m_grid->gridX);
  } else {
    temp = new OneDimCubicSplineFirst<T>(m_grid->gridX);
  }
  int nyz=nY*nZ;
  // Loop over all y and z
  for(int iy = 0; iy < nY; iy++){
    for(int iz = 0; iz < nZ; iz++){
      int first(index(0,iy,iz));
      T* restrict s = temp->data();
      for(int ix=0,cur=first; ix< nX; ix++,cur+=nyz,s++) *s = F[cur][source];

      temp->spline();

      const T* restrict sfit = temp->data(1);
      for(int ix=0,cur=first; ix< nX; ix++,cur+=nyz) F[cur][target] = *sfit++;
    }
  }

  delete temp;
}

template<class T>
void TriCubicSplineT<T>::UpdateY(int source, int target, bool periodic){

  typename OneDimGridFunctor<T>* temp=0;
  if(periodic) {
    temp = new OneDimCubicSplinePBC<T>(m_grid->gridY);
  } else {
    temp = new OneDimCubicSplineFirst<T>(m_grid->gridY);
  }

  // Loop over all x and z
  for(int ix = 0; ix < nX; ix++){
    for(int iz = 0; iz < nZ; iz++){
      int first(index(ix,0,iz));
      T* restrict s = temp->data();
      for(int iy=0,cur=first; iy< nY; iy++,cur+=nZ,s++) *s = F[cur][source];
      temp->spline();
      const T* restrict sfit = temp->data(1);
      for(int iy=0,cur=first; iy< nY; iy++,cur+=nZ) F[cur][target] = *sfit++;
    }
  }

  delete temp;
}

template<class T>
void TriCubicSplineT<T>::UpdateZ(int source, int target, bool periodic){

  typename OneDimGridFunctor<T>* temp=0;
  if(periodic) {
    temp = new OneDimCubicSplinePBC<T>(m_grid->gridZ);
  } else {
    temp = new OneDimCubicSplineFirst<T>(m_grid->gridZ);
  }

  // Loop over all x and z
  for(int ix = 0; ix < nX; ix++){
    for(int iy = 0; iy < nY; iy++){
      int first(index(ix,iy,0));
      T* restrict s = temp->data();
      for(int iz=0,cur=first; iz< nZ; iz++,cur++,s++) *s = F[cur][source];
      temp->spline();
      const T* restrict sfit = temp->data(1);
      for(int iz=0,cur=first; iz< nZ; iz++,cur++) F[cur][target] = *sfit++;
    }
  }

  delete temp;
}
#endif
