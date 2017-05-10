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
    
    



#ifndef QMCPLUSPLUS_XYZCUBICGRID_H
#define QMCPLUSPLUS_XYZCUBICGRID_H
#include <cmath>
#include "OhmmsPETE/TinyVector.h"
#include "Numerics/OneDimGridBase.h"

namespace qmcplusplus
{
/** Tri-cubic Splines with periodic boundary conditions and fixed first derivatives.
 *
 * Adapting TriCubicSpline implemented by K. Esler and D. Das.
 * Use stl containers
 */
template<typename T, typename Tg>
struct XYZCubicGrid
{

  typedef TinyVector<T,8>   KnotType;
  typedef OneDimGridBase<Tg> Grid1DType;

  /// functions which depend on the point where the interpolated value
  /// is required. t = (x - xi)/h
  inline Tg p1(Tg t)
  {
    return ((t-1.0)*(t-1.0)*(1.0+2.0*t));
  }
  inline Tg p2(Tg t)
  {
    return (t*t*(3.0-2.0*t));
  }
  inline Tg q1(Tg t)
  {
    return (t*(t-1.0)*(t-1.0));
  }
  inline Tg q2(Tg t)
  {
    return (t*t*(t-1.0));
  }
  inline Tg dp1(Tg t)
  {
    return (6.0*t*(t-1.0));
  }
  inline Tg dq1(Tg t)
  {
    return ((t-1.0)*(3.0*t-1.0));
  }
  inline Tg dp2(Tg t)
  {
    return (-dp1(t));
  }
  inline Tg dq2 (Tg t)
  {
    return ((3.0*t - 2.0)*t);
  }
  inline Tg d2p1(Tg t)
  {
    return (12.0*t-6.0);
  }
  inline Tg d2q1 (Tg t)
  {
    return (6.0*t - 4.0);
  }
  inline Tg d2p2 (Tg t)
  {
    return (-d2p1(t));
  }
  inline Tg d2q2 (Tg t)
  {
    return (6.0*t - 2.0);
  }

  bool OwnGrid;
  bool Periodic;
  int Loc;
  int ix, iy, iz;
  int nX, nY, nZ;
  Tg x_min, x_max, LengthX, OneOverLx;
  Tg y_min, y_max, LengthY, OneOverLy;
  Tg z_min, z_max, LengthZ, OneOverLz;
  Tg h,k,l,hinv,kinv,linv;
  Tg u,v,w;
  Tg a0,a1,a2,a3,b0,b1,b2,b3,c0,c1,c2,c3;
  Tg da0,da1,da2,da3,db0,db1,db2,db3,dc0,dc1,dc2,dc3;
  Tg d2a0,d2a1,d2a2,d2a3,d2b0,d2b1,d2b2,d2b3,d2c0,d2c1,d2c2,d2c3;

  T val, gradfX, gradfY, gradfZ, lapf;
  Grid1DType *gridX, *gridY, *gridZ;

  XYZCubicGrid(): OwnGrid(false),Loc(-1),gridX(0),gridY(0),gridZ(0) {}
  XYZCubicGrid(Grid1DType *xgrid, Grid1DType *ygrid,
               Grid1DType *zgrid)
  {
    setGridXYZ(xgrid,ygrid,zgrid);
  }

  inline
  void
  setGridXYZ(Grid1DType *xgrid, Grid1DType *ygrid, Grid1DType *zgrid)
  {
    gridX=xgrid;
    gridY=ygrid;
    gridZ=zgrid;
    x_min=gridX->rmin();
    x_max=gridX->rmax();
    LengthX=x_max-x_min;
    OneOverLx=1.0/LengthX;
    y_min=gridY->rmin();
    y_max=gridY->rmax();
    LengthY=y_max-y_min;
    OneOverLy=1.0/LengthY;
    z_min=gridZ->rmin();
    z_max=gridZ->rmax();
    LengthZ=z_max-z_min;
    OneOverLz=1.0/LengthZ;
    nX = gridX->size();
    nY = gridY->size();
    nZ = gridZ->size();
  }

  inline int size()
  {
    return nX*nY*nZ;
  }

  /** process xmlnode to set the grids
   */
  bool put(xmlNodePtr cur)
  {
    std::vector<Tg> ri(3,-5.0);
    std::vector<Tg> rf(3,5.0);
    std::vector<int> npts(3,101);
    cur = cur->xmlChildrenNode;
    int idir(0);
    while(cur != NULL)
    {
      std::string cname((const char*)(cur->name));
      if(cname == "grid")
      {
        const xmlChar* a=xmlGetProp(cur,(const xmlChar*)"dir");
        if(a)
        {
          idir=atoi((const char*)a);
        }
        a=xmlGetProp(cur,(const xmlChar*)"ri");
        if(a)
        {
          ri[idir]=atof((const char*)a);
        }
        a=xmlGetProp(cur,(const xmlChar*)"rf");
        if(a)
        {
          rf[idir]=atof((const char*)a);
        }
        a=xmlGetProp(cur,(const xmlChar*)"npts");
        if(a)
        {
          npts[idir]=atoi((const char*)a);
        }
      }
      cur=cur->next;
    }
    if(gridX ==0)
      gridX=new LinearGrid<Tg>;
    if(gridY ==0)
      gridY=new LinearGrid<Tg>;
    if(gridZ ==0)
      gridZ=new LinearGrid<Tg>;
    gridX->set(ri[0],rf[0],npts[0]);
    gridY->set(ri[1],rf[1],npts[1]);
    gridZ->set(ri[2],rf[2],npts[2]);
    setGridXYZ(gridX,gridY,gridZ);
    return true;
  }

  inline void setBC(bool pbc)
  {
    Periodic=pbc;
  }

  inline int index(int i, int j, int k) const
  {
    return k+nZ*(j+nY*i);
  }

  /** locate the grid point (x,y,z) and update the coefficients
   */
  inline void locate(Tg x, Tg y, Tg z, bool updateall)
  {
    //grid(X,Y,Z)->locate(r) evaluates the factors used by interpolations
    Loc=-1;
    if(Periodic)
    {
      x-=LengthX*std::floor(x*OneOverLx);
      y-=LengthY*std::floor(y*OneOverLy);
      z-=LengthZ*std::floor(z*OneOverLz);
      //T x0=std::fmod(x-x_min,LengthX);
      //x=x0-LengthX*static_cast<int>(x0*TwoOverLx)+x_min;
      //T y0=std::fmod(y-y_min,LengthY);
      //y=y0-LengthY*static_cast<int>(y0*TwoOverLy)+y_min;
      //T z0=std::fmod(z-z_min,LengthZ);
      //z=z0-LengthZ*static_cast<int>(z0*TwoOverLz)+z_min;
      //if(x<x_min) x+=LengthX;
      //else if(x>=x_max) x-=LengthX;
      //if(y<y_min) y+=LengthY;
      //else if(y>=y_max) y-=LengthY;
      //if(z<z_min) z+=LengthZ;
      //else if(z>=z_max) z-=LengthZ;
    }
    else
    {
      if(x<x_min || x > x_max)
        return;
      if(y<y_min || y > y_max)
        return;
      if(z<z_min || z > z_max)
        return;
    }
    ix=gridX->getIndex(x);
    iy=gridY->getIndex(y);
    iz=gridY->getIndex(z);
    h = gridX->dr(ix);
    hinv = 1.0/h;
    u = (x - gridX->r(ix))*hinv;
    k = gridY->dr(iy);
    kinv = 1.0/k;
    v = (y - gridY->r(iy))*kinv;
    l = gridZ->dr(iz);
    linv = 1.0/l;
    w = (z - gridZ->r(iz))*linv;
    Loc=index(ix,iy,iz);
    a0=p1(u);
    a1=p2(u);
    a2=h*q1(u);
    a3=h*q2(u);
    b0=p1(v);
    b1=p2(v);
    b2=k*q1(v);
    b3=k*q2(v);
    c0=p1(w);
    c1=p2(w);
    c2=l*q1(w);
    c3=l*q2(w);
    if(updateall)
    {
      da0=hinv*dp1(u);
      da1=hinv*dp2(u);
      da2=dq1(u);
      da3=dq2(u);
      db0=kinv*dp1(v);
      db1=kinv*dp2(v);
      db2=dq1(v);
      db3=dq2(v);
      dc0=linv*dp1(w);
      dc1=linv*dp2(w);
      dc2=dq1(w);
      dc3=dq2(w);
      d2a0=hinv*hinv*d2p1(u);
      d2a1=hinv*hinv*d2p2(u);
      d2a2=hinv*d2q1(u);
      d2a3=hinv*d2q2(u);
      d2b0=kinv*kinv*d2p1(v);
      d2b1=kinv*kinv*d2p2(v);
      d2b2=kinv*d2q1(v);
      d2b3=kinv*d2q2(v);
      d2c0=linv*linv*d2p1(w);
      d2c1=linv*linv*d2p2(w);
      d2c2=linv*d2q1(w);
      d2c3=linv*d2q2(w);
    }
  }

  //inline int update() {
  //  a0=p1(u); a1=p2(u); a2=h*q1(u); a3=h*q2(u);
  //  b0=p1(v); b1=p2(v); b2=k*q1(v); b3=k*q2(v);
  //  c0=p1(w); c1=p2(w); c2=l*q1(w); c3=l*q2(w);
  //  return Loc;
  //}

  inline void update(bool all)
  {
    if(Loc<0)
      return;
    a0=p1(u);
    a1=p2(u);
    a2=h*q1(u);
    a3=h*q2(u);
    b0=p1(v);
    b1=p2(v);
    b2=k*q1(v);
    b3=k*q2(v);
    c0=p1(w);
    c1=p2(w);
    c2=l*q1(w);
    c3=l*q2(w);
    if(all)
    {
      da0=hinv*dp1(u);
      da1=hinv*dp2(u);
      da2=dq1(u);
      da3=dq2(u);
      db0=kinv*dp1(v);
      db1=kinv*dp2(v);
      db2=dq1(v);
      db3=dq2(v);
      dc0=linv*dp1(w);
      dc1=linv*dp2(w);
      dc2=dq1(w);
      dc3=dq2(w);
      d2a0=hinv*hinv*d2p1(u);
      d2a1=hinv*hinv*d2p2(u);
      d2a2=hinv*d2q1(u);
      d2a3=hinv*d2q2(u);
      d2b0=kinv*kinv*d2p1(v);
      d2b1=kinv*kinv*d2p2(v);
      d2b2=kinv*d2q1(v);
      d2b3=kinv*d2q2(v);
      d2c0=linv*linv*d2p1(w);
      d2c1=linv*linv*d2p2(w);
      d2c2=linv*d2q1(w);
      d2c3=linv*d2q2(w);
    }
  }

  inline T evaluate(const KnotType& f000, const KnotType& f001,
                    const KnotType& f010, const KnotType& f011,
                    const KnotType& f100, const KnotType& f101,
                    const KnotType& f110, const KnotType& f111)
  {
    return
      a0*
      (b0*(f000[0]*c0+f001[0]*c1+f000[3]*c2+f001[3]*c3) +
       b1*(f010[0]*c0+f011[0]*c1+f010[3]*c2+f011[3]*c3) +
       b2*(f000[2]*c0+f001[2]*c1+f000[6]*c2+f001[6]*c3) +
       b3*(f010[2]*c0+f011[2]*c1+f010[6]*c2+f011[6]*c3))+
      a1 *
      (b0*(f100[0]*c0+f101[0]*c1+f100[3]*c2+f101[3]*c3) +
       b1*(f110[0]*c0+f111[0]*c1+f110[3]*c2+f111[3]*c3) +
       b2*(f100[2]*c0+f101[2]*c1+f100[6]*c2+f101[6]*c3) +
       b3*(f110[2]*c0+f111[2]*c1+f110[6]*c2+f111[6]*c3))+
      a2 *
      (b0*(f000[1]*c0+f001[1]*c1+f000[5]*c2+f001[5]*c3) +
       b1*(f010[1]*c0+f011[1]*c1+f010[5]*c2+f011[5]*c3) +
       b2*(f000[4]*c0+f001[4]*c1+f000[7]*c2+f001[7]*c3) +
       b3*(f010[4]*c0+f011[4]*c1+f010[7]*c2+f011[7]*c3))+
      a3 *
      (b0*(f100[1]*c0+f101[1]*c1+f100[5]*c2+f101[5]*c3) +
       b1*(f110[1]*c0+f111[1]*c1+f110[5]*c2+f111[5]*c3) +
       b2*(f100[4]*c0+f101[4]*c1+f100[7]*c2+f101[7]*c3) +
       b3*(f110[4]*c0+f111[4]*c1+f110[7]*c2+f111[7]*c3));
  }

  inline void evaluateAll(const KnotType& f000, const KnotType& f001,
                          const KnotType& f010, const KnotType& f011,
                          const KnotType& f100, const KnotType& f101,
                          const KnotType& f110, const KnotType& f111)
  {
    T Y000(f000[0]);    //   F
    T Y200(f000[1]);    //  dF/dx
    T Y020(f000[2]);    //  dF/dy
    T Y002(f000[3]);    //  dF/dz
    T Y220(f000[4]);    // d2F/dxdy
    T Y202(f000[5]);    // d2F/dxdz
    T Y022(f000[6]);    // d2F/dydz
    T Y222(f000[7]);    // d3F/dxdydz
    T Y001(f001[0]);    //   F
    T Y201(f001[1]);    //  dF/dx
    T Y021(f001[2]);    //  dF/dy
    T Y003(f001[3]);    //  dF/dz
    T Y221(f001[4]);    // d2F/dxdy
    T Y203(f001[5]);    // d2F/dxdz
    T Y023(f001[6]);    // d2F/dydz
    T Y223(f001[7]);    // d3F/dxdydz
    T Y010(f010[0]);    //   F
    T Y210(f010[1]);    //  dF/dx
    T Y030(f010[2]);    //  dF/dy
    T Y012(f010[3]);    //  dF/dz
    T Y230(f010[4]);    // d2F/dxdy
    T Y212(f010[5]);    // d2F/dxdz
    T Y032(f010[6]);    // d2F/dydz
    T Y232(f010[7]);    // d3F/dxdydz
    T Y011(f011[0]);  //   F
    T Y211(f011[1]);  //  dF/dx
    T Y031(f011[2]);  //  dF/dy
    T Y013(f011[3]);  //  dF/dz
    T Y231(f011[4]);  // d2F/dxdy
    T Y213(f011[5]);  // d2F/dxdz
    T Y033(f011[6]);  // d2F/dydz
    T Y233(f011[7]);  // d3F/dxdydz
    T Y100(f100[0]);      //   F
    T Y300(f100[1]);      //  dF/dx
    T Y120(f100[2]);      //  dF/dy
    T Y102(f100[3]);      //  dF/dz
    T Y320(f100[4]);      // d2F/dxdy
    T Y302(f100[5]);      // d2F/dxdz
    T Y122(f100[6]);      // d2F/dydz
    T Y322(f100[7]);      // d3F/dxdydz
    T Y101(f101[0]);    //   F
    T Y301(f101[1]);    //  dF/dx
    T Y121(f101[2]);    //  dF/dy
    T Y103(f101[3]);    //  dF/dz
    T Y321(f101[4]);    // d2F/dxdy
    T Y303(f101[5]);    // d2F/dxdz
    T Y123(f101[6]);    // d2F/dydz
    T Y323(f101[7]);    // d3F/dxdydz
    T Y110(f110[0]);    //   F
    T Y310(f110[1]);    //  dF/dx
    T Y130(f110[2]);    //  dF/dy
    T Y112(f110[3]);    //  dF/dz
    T Y330(f110[4]);    // d2F/dxdy
    T Y312(f110[5]);    // d2F/dxdz
    T Y132(f110[6]);    // d2F/dydz
    T Y332(f110[7]);    // d3F/dxdydz
    T Y111(f111[0]);  //   F
    T Y311(f111[1]);  //  dF/dx
    T Y131(f111[2]);  //  dF/dy
    T Y113(f111[3]);  //  dF/dz
    T Y331(f111[4]);  // d2F/dxdy
    T Y313(f111[5]);  // d2F/dxdz
    T Y133(f111[6]);  // d2F/dydz
    T Y333(f111[7]);  // d3F/dxdydz
    T term0 =
      b0*(Y000*c0+Y001*c1+Y002*c2+Y003*c3) +
      b1*(Y010*c0+Y011*c1+Y012*c2+Y013*c3) +
      b2*(Y020*c0+Y021*c1+Y022*c2+Y023*c3) +
      b3*(Y030*c0+Y031*c1+Y032*c2+Y033*c3);
    T term1 =
      b0*(Y100*c0+Y101*c1+Y102*c2+Y103*c3) +
      b1*(Y110*c0+Y111*c1+Y112*c2+Y113*c3) +
      b2*(Y120*c0+Y121*c1+Y122*c2+Y123*c3) +
      b3*(Y130*c0+Y131*c1+Y132*c2+Y133*c3);
    T term2 =
      b0*(Y200*c0+Y201*c1+Y202*c2+Y203*c3) +
      b1*(Y210*c0+Y211*c1+Y212*c2+Y213*c3) +
      b2*(Y220*c0+Y221*c1+Y222*c2+Y223*c3) +
      b3*(Y230*c0+Y231*c1+Y232*c2+Y233*c3);
    T term3 =
      b0*(Y300*c0+Y301*c1+Y302*c2+Y303*c3) +
      b1*(Y310*c0+Y311*c1+Y312*c2+Y313*c3) +
      b2*(Y320*c0+Y321*c1+Y322*c2+Y323*c3) +
      b3*(Y330*c0+Y331*c1+Y332*c2+Y333*c3);
    val = a0*term0 + a1*term1 + a2*term2 + a3*term3;
    gradfX = da0*term0 + da1*term1 + da2*term2 + da3*term3;
    gradfY =
      a0*
      (db0*(Y000*c0+Y001*c1+Y002*c2+Y003*c3) +
       db1*(Y010*c0+Y011*c1+Y012*c2+Y013*c3) +
       db2*(Y020*c0+Y021*c1+Y022*c2+Y023*c3) +
       db3*(Y030*c0+Y031*c1+Y032*c2+Y033*c3))+
      a1 *
      (db0*(Y100*c0+Y101*c1+Y102*c2+Y103*c3) +
       db1*(Y110*c0+Y111*c1+Y112*c2+Y113*c3) +
       db2*(Y120*c0+Y121*c1+Y122*c2+Y123*c3) +
       db3*(Y130*c0+Y131*c1+Y132*c2+Y133*c3))+
      a2 *
      (db0*(Y200*c0+Y201*c1+Y202*c2+Y203*c3) +
       db1*(Y210*c0+Y211*c1+Y212*c2+Y213*c3) +
       db2*(Y220*c0+Y221*c1+Y222*c2+Y223*c3) +
       db3*(Y230*c0+Y231*c1+Y232*c2+Y233*c3))+
      a3 *
      (db0*(Y300*c0+Y301*c1+Y302*c2+Y303*c3) +
       db1*(Y310*c0+Y311*c1+Y312*c2+Y313*c3) +
       db2*(Y320*c0+Y321*c1+Y322*c2+Y323*c3) +
       db3*(Y330*c0+Y331*c1+Y332*c2+Y333*c3));
    gradfZ =
      a0*
      (b0*(Y000*dc0+Y001*dc1+Y002*dc2+Y003*dc3) +
       b1*(Y010*dc0+Y011*dc1+Y012*dc2+Y013*dc3) +
       b2*(Y020*dc0+Y021*dc1+Y022*dc2+Y023*dc3) +
       b3*(Y030*dc0+Y031*dc1+Y032*dc2+Y033*dc3))+
      a1 *
      (b0*(Y100*dc0+Y101*dc1+Y102*dc2+Y103*dc3) +
       b1*(Y110*dc0+Y111*dc1+Y112*dc2+Y113*dc3) +
       b2*(Y120*dc0+Y121*dc1+Y122*dc2+Y123*dc3) +
       b3*(Y130*dc0+Y131*dc1+Y132*dc2+Y133*dc3))+
      a2 *
      (b0*(Y200*dc0+Y201*dc1+Y202*dc2+Y203*dc3) +
       b1*(Y210*dc0+Y211*dc1+Y212*dc2+Y213*dc3) +
       b2*(Y220*dc0+Y221*dc1+Y222*dc2+Y223*dc3) +
       b3*(Y230*dc0+Y231*dc1+Y232*dc2+Y233*dc3))+
      a3 *
      (b0*(Y300*dc0+Y301*dc1+Y302*dc2+Y303*dc3) +
       b1*(Y310*dc0+Y311*dc1+Y312*dc2+Y313*dc3) +
       b2*(Y320*dc0+Y321*dc1+Y322*dc2+Y323*dc3) +
       b3*(Y330*dc0+Y331*dc1+Y332*dc2+Y333*dc3));
    lapf =
      d2a0*
      (b0*(Y000*c0+Y001*c1+Y002*c2+Y003*c3) +
       b1*(Y010*c0+Y011*c1+Y012*c2+Y013*c3) +
       b2*(Y020*c0+Y021*c1+Y022*c2+Y023*c3) +
       b3*(Y030*c0+Y031*c1+Y032*c2+Y033*c3))+
      d2a1 *
      (b0*(Y100*c0+Y101*c1+Y102*c2+Y103*c3) +
       b1*(Y110*c0+Y111*c1+Y112*c2+Y113*c3) +
       b2*(Y120*c0+Y121*c1+Y122*c2+Y123*c3) +
       b3*(Y130*c0+Y131*c1+Y132*c2+Y133*c3))+
      d2a2 *
      (b0*(Y200*c0+Y201*c1+Y202*c2+Y203*c3) +
       b1*(Y210*c0+Y211*c1+Y212*c2+Y213*c3) +
       b2*(Y220*c0+Y221*c1+Y222*c2+Y223*c3) +
       b3*(Y230*c0+Y231*c1+Y232*c2+Y233*c3))+
      d2a3 *
      (b0*(Y300*c0+Y301*c1+Y302*c2+Y303*c3) +
       b1*(Y310*c0+Y311*c1+Y312*c2+Y313*c3) +
       b2*(Y320*c0+Y321*c1+Y322*c2+Y323*c3) +
       b3*(Y330*c0+Y331*c1+Y332*c2+Y333*c3))+
      a0*
      (d2b0*(Y000*c0+Y001*c1+Y002*c2+Y003*c3) +
       d2b1*(Y010*c0+Y011*c1+Y012*c2+Y013*c3) +
       d2b2*(Y020*c0+Y021*c1+Y022*c2+Y023*c3) +
       d2b3*(Y030*c0+Y031*c1+Y032*c2+Y033*c3))+
      a1 *
      (d2b0*(Y100*c0+Y101*c1+Y102*c2+Y103*c3) +
       d2b1*(Y110*c0+Y111*c1+Y112*c2+Y113*c3) +
       d2b2*(Y120*c0+Y121*c1+Y122*c2+Y123*c3) +
       d2b3*(Y130*c0+Y131*c1+Y132*c2+Y133*c3))+
      a2 *
      (d2b0*(Y200*c0+Y201*c1+Y202*c2+Y203*c3) +
       d2b1*(Y210*c0+Y211*c1+Y212*c2+Y213*c3) +
       d2b2*(Y220*c0+Y221*c1+Y222*c2+Y223*c3) +
       d2b3*(Y230*c0+Y231*c1+Y232*c2+Y233*c3))+
      a3 *
      (d2b0*(Y300*c0+Y301*c1+Y302*c2+Y303*c3) +
       d2b1*(Y310*c0+Y311*c1+Y312*c2+Y313*c3) +
       d2b2*(Y320*c0+Y321*c1+Y322*c2+Y323*c3) +
       d2b3*(Y330*c0+Y331*c1+Y332*c2+Y333*c3)) +
      a0*
      (b0*(Y000*d2c0+Y001*d2c1+Y002*d2c2+Y003*d2c3) +
       b1*(Y010*d2c0+Y011*d2c1+Y012*d2c2+Y013*d2c3) +
       b2*(Y020*d2c0+Y021*d2c1+Y022*d2c2+Y023*d2c3) +
       b3*(Y030*d2c0+Y031*d2c1+Y032*d2c2+Y033*d2c3))+
      a1 *
      (b0*(Y100*d2c0+Y101*d2c1+Y102*d2c2+Y103*d2c3) +
       b1*(Y110*d2c0+Y111*d2c1+Y112*d2c2+Y113*d2c3) +
       b2*(Y120*d2c0+Y121*d2c1+Y122*d2c2+Y123*d2c3) +
       b3*(Y130*d2c0+Y131*d2c1+Y132*d2c2+Y133*d2c3))+
      a2 *
      (b0*(Y200*d2c0+Y201*d2c1+Y202*d2c2+Y203*d2c3) +
       b1*(Y210*d2c0+Y211*d2c1+Y212*d2c2+Y213*d2c3) +
       b2*(Y220*d2c0+Y221*d2c1+Y222*d2c2+Y223*d2c3) +
       b3*(Y230*d2c0+Y231*d2c1+Y232*d2c2+Y233*d2c3))+
      a3 *
      (b0*(Y300*d2c0+Y301*d2c1+Y302*d2c2+Y303*d2c3) +
       b1*(Y310*d2c0+Y311*d2c1+Y312*d2c2+Y313*d2c3) +
       b2*(Y320*d2c0+Y321*d2c1+Y322*d2c2+Y323*d2c3) +
       b3*(Y330*d2c0+Y331*d2c1+Y332*d2c2+Y333*d2c3));
  }
};
}
#endif
