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
    
    




template<typename T, bool PBC0, bool PBC1, bool PBC2>
inline TricubicBsplineGridBC<T,PBC0,PBC1,PBC2>::TricubicBsplineGridBC()
{
  A(0,0) = -1.0/6.0;
  A(0,1) =  3.0/6.0;
  A(0,2) = -3.0/6.0;
  A(0,3) = 1.0/6.0;
  A(1,0) =  3.0/6.0;
  A(1,1) = -6.0/6.0;
  A(1,2) =  3.0/6.0;
  A(1,3) = 0.0/6.0;
  A(2,0) = -3.0/6.0;
  A(2,1) =  0.0/6.0;
  A(2,2) =  3.0/6.0;
  A(2,3) = 0.0/6.0;
  A(3,0) =  1.0/6.0;
  A(3,1) =  4.0/6.0;
  A(3,2) =  1.0/6.0;
  A(3,3) = 0.0/6.0;
  px[3] = 1.0;
  py[3] = 1.0;
  pz[3] = 1.0;
  dA(0,0)= 0.0;
  dA(0,1)= 0.0;
  dA(0,2)= 0.0;
  dA(0,3)= 0.0;
  dA(1,0)=-0.5;
  dA(1,1)= 1.5;
  dA(1,2)=-1.5;
  dA(1,3)= 0.5;
  dA(2,0)= 1.0;
  dA(2,1)=-2.0;
  dA(2,2)= 1.0;
  dA(2,3)= 0.0;
  dA(3,0)=-0.5;
  dA(3,1)= 0.0;
  dA(3,2)= 0.5;
  dA(3,3)= 0.0;
  d2A(0,0)= 0.0;
  d2A(0,1)= 0.0;
  d2A(0,2)= 0.0;
  d2A(0,3)= 0.0;
  d2A(1,0)= 0.0;
  d2A(1,1)= 0.0;
  d2A(1,2)= 0.0;
  d2A(1,3)= 0.0;
  d2A(2,0)=-1.0;
  d2A(2,1)= 3.0;
  d2A(2,2)=-3.0;
  d2A(2,3)= 1.0;
  d2A(3,0)= 1.0;
  d2A(3,1)=-2.0;
  d2A(3,2)= 1.0;
  d2A(3,3)= 0.0;
  d3A(0,0)= 0.0;
  d3A(0,1)= 0.0;
  d3A(0,2)= 0.0;
  d3A(0,3)= 0.0;
  d3A(1,0)= 0.0;
  d3A(1,1)= 0.0;
  d3A(1,2)= 0.0;
  d3A(1,3)= 0.0;
  d3A(2,0)= 0.0;
  d3A(2,1)= 0.0;
  d3A(1,2)= 2.0;
  d3A(2,3)= 0.0;
  d3A(3,0)=-1.0;
  d3A(3,1)= 3.0;
  d3A(3,2)=-3.0;
  d3A(3,3)= 1.0;
}

template<typename T, bool PBC0, bool PBC1, bool PBC2>
inline TricubicBsplineGridBC<T,PBC0,PBC1,PBC2>::TricubicBsplineGridBC(const TricubicBsplineGridBC<T,PBC0,PBC1,PBC2>& rhs)
{
  setGrid(rhs.BC0.Lower,rhs.BC0.Upper,
          rhs.BC1.Lower, rhs.BC1.Upper,
          rhs.BC2.Lower, rhs.BC2.Upper,
          rhs.Nx, rhs.Ny, rhs.Nz, true);
}

template<typename T, bool PBC0, bool PBC1, bool PBC2>
inline TricubicBsplineGridBC<T,PBC0,PBC1,PBC2>&
TricubicBsplineGridBC<T,PBC0,PBC1,PBC2>::operator=(const TricubicBsplineGridBC<T,PBC0,PBC1,PBC2>& rhs)
{
  if(this != &rhs)//prevent self-assignment
  {
    setGrid(rhs.BC0.Lower,rhs.BC0.Upper,
            rhs.BC1.Lower, rhs.BC1.Upper,
            rhs.BC2.Lower, rhs.BC2.Upper,
            rhs.Nx, rhs.Ny, rhs.Nz, true);
  }
  return *this;
}

template<typename T, bool PBC0, bool PBC1, bool PBC2>
inline void TricubicBsplineGridBC<T,PBC0,PBC1,PBC2>::setGrid(real_type xi, real_type xf,
    real_type yi, real_type yf, real_type zi, real_type zf,
    int nx, int ny, int nz,
    bool openend)
{
  Interpolating=true;
  Periodic=(PBC0||PBC1||PBC2);//periodic if any of the three dimensions is periodic
  if(openend)
  {
    Nx=nx;
    Ny=ny;
    Nz=nz;
  }
  else
  {
    Nx=nx-1;
    Ny=ny-1;
    Nz=nz-1;
  }
  BC0.init(Nx,xi,xf);
  dx=BC0.Delta;
  dxInv=BC0.OneOverDelta;
  BC1.init(Ny,yi,yf);
  dy=BC1.Delta;
  dyInv=BC1.OneOverDelta;
  BC2.init(Nz,zi,zf);
  dz=BC2.Delta;
  dzInv=BC2.OneOverDelta;
}

template<typename T, bool PBC0, bool PBC1, bool PBC2>
inline bool TricubicBsplineGridBC<T,PBC0,PBC1,PBC2>::Find(real_type x, real_type y, real_type z)
{
  if(BC0.outofbound(x,ix0))
    return false;
  ix1=ix0+1;
  ix2=ix0+2;
  ix3=ix0+3;
  if(BC1.outofbound(y,iy0))
    return false;
  iy1=iy0+1;
  iy2=iy0+2;
  iy3=iy0+3;
  if(BC2.outofbound(z,iz0))
    return false;
  iz1=iz0+1;
  iz2=iz0+2;
  iz3=iz0+3;
  px[0] = x*x*x;
  py[0] = y*y*y;
  pz[0] = z*z*z;
  px[1] = x*x;
  py[1] = y*y;
  pz[1] = z*z;
  px[2] = x;
  py[2] = y;
  pz[2] = z;
  //px[3] = 1.0;    py[3] = 1.0;   pz[3] = 1.0;
  a=dot(px,A);
  b=dot(py,A);
  c=dot(pz,A);
  return true;
}

template<typename T, bool PBC0, bool PBC1, bool PBC2>
inline bool TricubicBsplineGridBC<T,PBC0,PBC1,PBC2>::FindAll(real_type x, real_type y, real_type z)
{
  if(Find(x,y,z))
  {
    //da=dot(px,dA);
    //db=dot(py,dA);
    //dc=dot(pz,dA);
    //// First derivatives;
    da[0] = dA(1,0)*px[1];
    da[1] =dA(1,1)*px[1];
    da[2] =dA(1,2)*px[1];
    da[3] =dA(1,3)*px[1];
    da[0]+= dA(2,0)*px[2];
    da[1]+=dA(2,1)*px[2];
    da[2]+=dA(2,2)*px[2];
    da[3]+=dA(2,3)*px[2];
    da[0]+= dA(3,0)*px[3];
    da[1]+=dA(3,1)*px[3];
    da[2]+=dA(3,2)*px[3];
    da[3]+=dA(3,3)*px[3];
    db[0] = dA(1,0)*py[1];
    db[1] =dA(1,1)*py[1];
    db[2] =dA(1,2)*py[1];
    db[3] =dA(1,3)*py[1];
    db[0]+= dA(2,0)*py[2];
    db[1]+=dA(2,1)*py[2];
    db[2]+=dA(2,2)*py[2];
    db[3]+=dA(2,3)*py[2];
    db[0]+= dA(3,0)*py[3];
    db[1]+=dA(3,1)*py[3];
    db[2]+=dA(3,2)*py[3];
    db[3]+=dA(3,3)*py[3];
    dc[0] = dA(1,0)*pz[1];
    dc[1] =dA(1,1)*pz[1];
    dc[2] =dA(1,2)*pz[1];
    dc[3] =dA(1,3)*pz[1];
    dc[0]+= dA(2,0)*pz[2];
    dc[1]+=dA(2,1)*pz[2];
    dc[2]+=dA(2,2)*pz[2];
    dc[3]+=dA(2,3)*pz[2];
    dc[0]+= dA(3,0)*pz[3];
    dc[1]+=dA(3,1)*pz[3];
    dc[2]+=dA(3,2)*pz[3];
    dc[3]+=dA(3,3)*pz[3];
    // Second derivatives
    //d2a=dot(px,d2A);
    //d2b=dot(py,d2A);
    //d2c=dot(pz,d2A);
    d2a[0] = d2A(2,0)*px[2];
    d2a[1] =d2A(2,1)*px[2];
    d2a[2] =d2A(2,2)*px[2];
    d2a[3] =d2A(2,3)*px[2];
    d2a[0]+= d2A(3,0)*px[3];
    d2a[1]+=d2A(3,1)*px[3];
    d2a[2]+=d2A(3,2)*px[3];
    d2a[3]+=d2A(3,3)*px[3];
    d2b[0] = d2A(2,0)*py[2];
    d2b[1] =d2A(2,1)*py[2];
    d2b[2] =d2A(2,2)*py[2];
    d2b[3] =d2A(2,3)*py[2];
    d2b[0]+= d2A(3,0)*py[3];
    d2b[1]+=d2A(3,1)*py[3];
    d2b[2]+=d2A(3,2)*py[3];
    d2b[3]+=d2A(3,3)*py[3];
    d2c[0] = d2A(2,0)*pz[2];
    d2c[1] =d2A(2,1)*pz[2];
    d2c[2] =d2A(2,2)*pz[2];
    d2c[3] =d2A(2,3)*pz[2];
    d2c[0]+= d2A(3,0)*pz[3];
    d2c[1]+=d2A(3,1)*pz[3];
    d2c[2]+=d2A(3,2)*pz[3];
    d2c[3]+=d2A(3,3)*pz[3];
    return true;
  }
  else
    return false;
}

template<typename T, bool PBC0, bool PBC1, bool PBC2>
inline T TricubicBsplineGridBC<T,PBC0,PBC1,PBC2>::evaluate(const Array<T,3>& P) const
{
  return
    (a[0]*(b[0]*(c[0]*P(ix0,iy0,iz0)+c[1]*P(ix0,iy0,iz1)+c[2]*P(ix0,iy0,iz2)+c[3]*P(ix0,iy0,iz3))+
           b[1]*(c[0]*P(ix0,iy1,iz0)+c[1]*P(ix0,iy1,iz1)+c[2]*P(ix0,iy1,iz2)+c[3]*P(ix0,iy1,iz3))+
           b[2]*(c[0]*P(ix0,iy2,iz0)+c[1]*P(ix0,iy2,iz1)+c[2]*P(ix0,iy2,iz2)+c[3]*P(ix0,iy2,iz3))+
           b[3]*(c[0]*P(ix0,iy3,iz0)+c[1]*P(ix0,iy3,iz1)+c[2]*P(ix0,iy3,iz2)+c[3]*P(ix0,iy3,iz3)))+
     a[1]*(b[0]*(c[0]*P(ix1,iy0,iz0)+c[1]*P(ix1,iy0,iz1)+c[2]*P(ix1,iy0,iz2)+c[3]*P(ix1,iy0,iz3))+
           b[1]*(c[0]*P(ix1,iy1,iz0)+c[1]*P(ix1,iy1,iz1)+c[2]*P(ix1,iy1,iz2)+c[3]*P(ix1,iy1,iz3))+
           b[2]*(c[0]*P(ix1,iy2,iz0)+c[1]*P(ix1,iy2,iz1)+c[2]*P(ix1,iy2,iz2)+c[3]*P(ix1,iy2,iz3))+
           b[3]*(c[0]*P(ix1,iy3,iz0)+c[1]*P(ix1,iy3,iz1)+c[2]*P(ix1,iy3,iz2)+c[3]*P(ix1,iy3,iz3)))+
     a[2]*(b[0]*(c[0]*P(ix2,iy0,iz0)+c[1]*P(ix2,iy0,iz1)+c[2]*P(ix2,iy0,iz2)+c[3]*P(ix2,iy0,iz3))+
           b[1]*(c[0]*P(ix2,iy1,iz0)+c[1]*P(ix2,iy1,iz1)+c[2]*P(ix2,iy1,iz2)+c[3]*P(ix2,iy1,iz3))+
           b[2]*(c[0]*P(ix2,iy2,iz0)+c[1]*P(ix2,iy2,iz1)+c[2]*P(ix2,iy2,iz2)+c[3]*P(ix2,iy2,iz3))+
           b[3]*(c[0]*P(ix2,iy3,iz0)+c[1]*P(ix2,iy3,iz1)+c[2]*P(ix2,iy3,iz2)+c[3]*P(ix2,iy3,iz3)))+
     a[3]*(b[0]*(c[0]*P(ix3,iy0,iz0)+c[1]*P(ix3,iy0,iz1)+c[2]*P(ix3,iy0,iz2)+c[3]*P(ix3,iy0,iz3))+
           b[1]*(c[0]*P(ix3,iy1,iz0)+c[1]*P(ix3,iy1,iz1)+c[2]*P(ix3,iy1,iz2)+c[3]*P(ix3,iy1,iz3))+
           b[2]*(c[0]*P(ix3,iy2,iz0)+c[1]*P(ix3,iy2,iz1)+c[2]*P(ix3,iy2,iz2)+c[3]*P(ix3,iy2,iz3))+
           b[3]*(c[0]*P(ix3,iy3,iz0)+c[1]*P(ix3,iy3,iz1)+c[2]*P(ix3,iy3,iz2)+c[3]*P(ix3,iy3,iz3))));
}

template<typename T, bool PBC0, bool PBC1, bool PBC2>
inline T TricubicBsplineGridBC<T,PBC0,PBC1,PBC2>::evaluate(const Array<T,3>& P,
    TinyVector<T,3>& grad, T& laplacian) const
{
  //put the blocks in cache
  TinyVector<T,4> p00(P(ix0,iy0,iz0),P(ix0,iy0,iz1),P(ix0,iy0,iz2),P(ix0,iy0,iz3));
  TinyVector<T,4> p01(P(ix0,iy1,iz0),P(ix0,iy1,iz1),P(ix0,iy1,iz2),P(ix0,iy1,iz3));
  TinyVector<T,4> p02(P(ix0,iy2,iz0),P(ix0,iy2,iz1),P(ix0,iy2,iz2),P(ix0,iy2,iz3));
  TinyVector<T,4> p03(P(ix0,iy3,iz0),P(ix0,iy3,iz1),P(ix0,iy3,iz2),P(ix0,iy3,iz3));
  TinyVector<T,4> p10(P(ix1,iy0,iz0),P(ix1,iy0,iz1),P(ix1,iy0,iz2),P(ix1,iy0,iz3));
  TinyVector<T,4> p11(P(ix1,iy1,iz0),P(ix1,iy1,iz1),P(ix1,iy1,iz2),P(ix1,iy1,iz3));
  TinyVector<T,4> p12(P(ix1,iy2,iz0),P(ix1,iy2,iz1),P(ix1,iy2,iz2),P(ix1,iy2,iz3));
  TinyVector<T,4> p13(P(ix1,iy3,iz0),P(ix1,iy3,iz1),P(ix1,iy3,iz2),P(ix1,iy3,iz3));
  TinyVector<T,4> p20(P(ix2,iy0,iz0),P(ix2,iy0,iz1),P(ix2,iy0,iz2),P(ix2,iy0,iz3));
  TinyVector<T,4> p21(P(ix2,iy1,iz0),P(ix2,iy1,iz1),P(ix2,iy1,iz2),P(ix2,iy1,iz3));
  TinyVector<T,4> p22(P(ix2,iy2,iz0),P(ix2,iy2,iz1),P(ix2,iy2,iz2),P(ix2,iy2,iz3));
  TinyVector<T,4> p23(P(ix2,iy3,iz0),P(ix2,iy3,iz1),P(ix2,iy3,iz2),P(ix2,iy3,iz3));
  TinyVector<T,4> p30(P(ix3,iy0,iz0),P(ix3,iy0,iz1),P(ix3,iy0,iz2),P(ix3,iy0,iz3));
  TinyVector<T,4> p31(P(ix3,iy1,iz0),P(ix3,iy1,iz1),P(ix3,iy1,iz2),P(ix3,iy1,iz3));
  TinyVector<T,4> p32(P(ix3,iy2,iz0),P(ix3,iy2,iz1),P(ix3,iy2,iz2),P(ix3,iy2,iz3));
  TinyVector<T,4> p33(P(ix3,iy3,iz0),P(ix3,iy3,iz1),P(ix3,iy3,iz2),P(ix3,iy3,iz3));
  // Save some operations by factorizing computation.
  Tensor<T,4> cP;
  //cP(0,0) = dot(c, p00);
  //cP(0,1) = dot(c, p01);
  //cP(0,2) = dot(c, p02);
  //cP(0,3) = dot(c, p03);
  //cP(1,0) = dot(c, p10);
  //cP(1,1) = dot(c, p11);
  //cP(1,2) = dot(c, p12);
  //cP(1,3) = dot(c, p13);
  //cP(2,0) = dot(c, p20);
  //cP(2,1) = dot(c, p21);
  //cP(2,2) = dot(c, p22);
  //cP(2,3) = dot(c, p23);
  //cP(3,0) = dot(c, p30);
  //cP(3,1) = dot(c, p31);
  //cP(3,2) = dot(c, p32);
  //cP(3,3) = dot(c, p33);
  cP(0,0) = dot(p00, c);
  cP(0,1) = dot(p01, c);
  cP(0,2) = dot(p02, c);
  cP(0,3) = dot(p03, c);
  cP(1,0) = dot(p10, c);
  cP(1,1) = dot(p11, c);
  cP(1,2) = dot(p12, c);
  cP(1,3) = dot(p13, c);
  cP(2,0) = dot(p20, c);
  cP(2,1) = dot(p21, c);
  cP(2,2) = dot(p22, c);
  cP(2,3) = dot(p23, c);
  cP(3,0) = dot(p30, c);
  cP(3,1) = dot(p31, c);
  cP(3,2) = dot(p32, c);
  cP(3,3) = dot(p33, c);
  TinyVector<T,4> bcP=dot(cP,b);
  // Compute value
  //T val = dot(a,bcP);
  // Compute gradient
  grad[0] = dxInv *dot(da,bcP);
  grad[1] = dyInv *
            (a[0]*(cP(0,0)*db[0]+cP(0,1)*db[1]+cP(0,2)*db[2]+cP(0,3)*db[3]) +
             a[1]*(cP(1,0)*db[0]+cP(1,1)*db[1]+cP(1,2)*db[2]+cP(1,3)*db[3]) +
             a[2]*(cP(2,0)*db[0]+cP(2,1)*db[1]+cP(2,2)*db[2]+cP(2,3)*db[3]) +
             a[3]*(cP(3,0)*db[0]+cP(3,1)*db[1]+cP(3,2)*db[2]+cP(3,3)*db[3]));
  grad[2] = dzInv *
            (a[0]*(b[0]*dot(dc, p00)+ b[1]*dot(dc, p01)+ b[2]*dot(dc, p02)+ b[3]*dot(dc, p03))+
             a[1]*(b[0]*dot(dc, p10)+ b[1]*dot(dc, p11)+ b[2]*dot(dc, p12)+ b[3]*dot(dc, p13))+
             a[2]*(b[0]*dot(dc, p20)+ b[1]*dot(dc, p21)+ b[2]*dot(dc, p22)+ b[3]*dot(dc, p23))+
             a[3]*(b[0]*dot(dc, p30)+ b[1]*dot(dc, p31)+ b[2]*dot(dc, p32)+ b[3]*dot(dc, p33)));
  // Compute laplacian
  laplacian = dxInv * dxInv *
              (d2a[0]*bcP[0] + d2a[1]*bcP[1] + d2a[2]*bcP[2] + d2a[3]*bcP[3]);
  laplacian += dyInv * dyInv *
               (a[0]*(cP(0,0)*d2b[0]+cP(0,1)*d2b[1]+cP(0,2)*d2b[2]+cP(0,3)*d2b[3]) +
                a[1]*(cP(1,0)*d2b[0]+cP(1,1)*d2b[1]+cP(1,2)*d2b[2]+cP(1,3)*d2b[3]) +
                a[2]*(cP(2,0)*d2b[0]+cP(2,1)*d2b[1]+cP(2,2)*d2b[2]+cP(2,3)*d2b[3]) +
                a[3]*(cP(3,0)*d2b[0]+cP(3,1)*d2b[1]+cP(3,2)*d2b[2]+cP(3,3)*d2b[3]));
  laplacian += dzInv * dzInv *
               (a[0]*(b[0]*dot(d2c, p00)+ b[1]*dot(d2c, p01)+ b[2]*dot(d2c, p02)+ b[3]*dot(d2c, p03))+
                a[1]*(b[0]*dot(d2c, p10)+ b[1]*dot(d2c, p11)+ b[2]*dot(d2c, p12)+ b[3]*dot(d2c, p13))+
                a[2]*(b[0]*dot(d2c, p20)+ b[1]*dot(d2c, p21)+ b[2]*dot(d2c, p22)+ b[3]*dot(d2c, p23))+
                a[3]*(b[0]*dot(d2c, p30)+ b[1]*dot(d2c, p31)+ b[2]*dot(d2c, p32)+ b[3]*dot(d2c, p33)));
  return dot(a,bcP);
}

template<typename T, bool PBC0, bool PBC1, bool PBC2>
inline T TricubicBsplineGridBC<T,PBC0,PBC1,PBC2>::evaluate(const Array<T,3>& P,
    TinyVector<T,3>& grad, Tensor<T,3>& secDerivs) const
{
  //put the blocks in cache
  TinyVector<T,4> p00(P(ix0,iy0,iz0),P(ix0,iy0,iz1),P(ix0,iy0,iz2),P(ix0,iy0,iz3));
  TinyVector<T,4> p01(P(ix0,iy1,iz0),P(ix0,iy1,iz1),P(ix0,iy1,iz2),P(ix0,iy1,iz3));
  TinyVector<T,4> p02(P(ix0,iy2,iz0),P(ix0,iy2,iz1),P(ix0,iy2,iz2),P(ix0,iy2,iz3));
  TinyVector<T,4> p03(P(ix0,iy3,iz0),P(ix0,iy3,iz1),P(ix0,iy3,iz2),P(ix0,iy3,iz3));
  TinyVector<T,4> p10(P(ix1,iy0,iz0),P(ix1,iy0,iz1),P(ix1,iy0,iz2),P(ix1,iy0,iz3));
  TinyVector<T,4> p11(P(ix1,iy1,iz0),P(ix1,iy1,iz1),P(ix1,iy1,iz2),P(ix1,iy1,iz3));
  TinyVector<T,4> p12(P(ix1,iy2,iz0),P(ix1,iy2,iz1),P(ix1,iy2,iz2),P(ix1,iy2,iz3));
  TinyVector<T,4> p13(P(ix1,iy3,iz0),P(ix1,iy3,iz1),P(ix1,iy3,iz2),P(ix1,iy3,iz3));
  TinyVector<T,4> p20(P(ix2,iy0,iz0),P(ix2,iy0,iz1),P(ix2,iy0,iz2),P(ix2,iy0,iz3));
  TinyVector<T,4> p21(P(ix2,iy1,iz0),P(ix2,iy1,iz1),P(ix2,iy1,iz2),P(ix2,iy1,iz3));
  TinyVector<T,4> p22(P(ix2,iy2,iz0),P(ix2,iy2,iz1),P(ix2,iy2,iz2),P(ix2,iy2,iz3));
  TinyVector<T,4> p23(P(ix2,iy3,iz0),P(ix2,iy3,iz1),P(ix2,iy3,iz2),P(ix2,iy3,iz3));
  TinyVector<T,4> p30(P(ix3,iy0,iz0),P(ix3,iy0,iz1),P(ix3,iy0,iz2),P(ix3,iy0,iz3));
  TinyVector<T,4> p31(P(ix3,iy1,iz0),P(ix3,iy1,iz1),P(ix3,iy1,iz2),P(ix3,iy1,iz3));
  TinyVector<T,4> p32(P(ix3,iy2,iz0),P(ix3,iy2,iz1),P(ix3,iy2,iz2),P(ix3,iy2,iz3));
  TinyVector<T,4> p33(P(ix3,iy3,iz0),P(ix3,iy3,iz1),P(ix3,iy3,iz2),P(ix3,iy3,iz3));
  // Save some operations by factorizing computation.
  Tensor<T,4> cP, dcP;
  cP(0,0) = dot(c, p00);
  cP(0,1) = dot(c, p01);
  cP(0,2) = dot(c, p02);
  cP(0,3) = dot(c, p03);
  cP(1,0) = dot(c, p10);
  cP(1,1) = dot(c, p11);
  cP(1,2) = dot(c, p12);
  cP(1,3) = dot(c, p13);
  cP(2,0) = dot(c, p20);
  cP(2,1) = dot(c, p21);
  cP(2,2) = dot(c, p22);
  cP(2,3) = dot(c, p23);
  cP(3,0) = dot(c, p30);
  cP(3,1) = dot(c, p31);
  cP(3,2) = dot(c, p32);
  cP(3,3) = dot(c, p33);
  dcP(0,0) = dot(dc, p00);
  dcP(0,1) = dot(dc, p01);
  dcP(0,2) = dot(dc, p02);
  dcP(0,3) = dot(dc, p03);
  dcP(1,0) = dot(dc, p10);
  dcP(1,1) = dot(dc, p11);
  dcP(1,2) = dot(dc, p12);
  dcP(1,3) = dot(dc, p13);
  dcP(2,0) = dot(dc, p20);
  dcP(2,1) = dot(dc, p21);
  dcP(2,2) = dot(dc, p22);
  dcP(2,3) = dot(dc, p23);
  dcP(3,0) = dot(dc, p30);
  dcP(3,1) = dot(dc, p31);
  dcP(3,2) = dot(dc, p32);
  dcP(3,3) = dot(dc, p33);
  TinyVector<T,4> bcP=dot(cP,b);
  TinyVector<T,4> bdcP=dot(dcP,b);
  // Compute value
  //val = a[0]*bcP[0] + a[1]*bcP[1] + a[2]*bcP[2] + a[3]*bcP[3];
  // Compute gradient
  grad[0] = dxInv * dot(da,bcP);
  grad[1] = dyInv *
            (a[0]*(cP(0,0)*db[0]+cP(0,1)*db[1]+cP(0,2)*db[2]+cP(0,3)*db[3]) +
             a[1]*(cP(1,0)*db[0]+cP(1,1)*db[1]+cP(1,2)*db[2]+cP(1,3)*db[3]) +
             a[2]*(cP(2,0)*db[0]+cP(2,1)*db[1]+cP(2,2)*db[2]+cP(2,3)*db[3]) +
             a[3]*(cP(3,0)*db[0]+cP(3,1)*db[1]+cP(3,2)*db[2]+cP(3,3)*db[3]));
  grad[2] = dzInv *  dot(a,bdcP);
  // Compute laplacian
  secDerivs(0,0) = dxInv * dxInv * dot(d2a,bcP);
  secDerivs(0,1) = secDerivs(1,0) = dxInv * dyInv *
                                    (da[0]*(cP(0,0)*db[0]+cP(0,1)*db[1]+cP(0,2)*db[2]+cP(0,3)*db[3]) +
                                     da[1]*(cP(1,0)*db[0]+cP(1,1)*db[1]+cP(1,2)*db[2]+cP(1,3)*db[3]) +
                                     da[2]*(cP(2,0)*db[0]+cP(2,1)*db[1]+cP(2,2)*db[2]+cP(2,3)*db[3]) +
                                     da[3]*(cP(3,0)*db[0]+cP(3,1)*db[1]+cP(3,2)*db[2]+cP(3,3)*db[3]));
  secDerivs(0,2) = secDerivs(2,0) = dxInv * dzInv *
                                    (da[0]*bdcP[0] + da[1]*bdcP[1] + da[2]*bdcP[2] + da[3]*bdcP[3]);
  secDerivs(1,1) = dyInv * dyInv *
                   (a[0]*(cP(0,0)*d2b[0]+cP(0,1)*d2b[1]+cP(0,2)*d2b[2]+cP(0,3)*d2b[3]) +
                    a[1]*(cP(1,0)*d2b[0]+cP(1,1)*d2b[1]+cP(1,2)*d2b[2]+cP(1,3)*d2b[3]) +
                    a[2]*(cP(2,0)*d2b[0]+cP(2,1)*d2b[1]+cP(2,2)*d2b[2]+cP(2,3)*d2b[3]) +
                    a[3]*(cP(3,0)*d2b[0]+cP(3,1)*d2b[1]+cP(3,2)*d2b[2]+cP(3,3)*d2b[3]));
  secDerivs(1,2) = secDerivs(2,1) = dyInv *dzInv *
                                    (a[0]*(db[0]*dcP(0,0) + db[1]*dcP(0,1) + db[2]*dcP(0,2) + db[3]*dcP(0,3))+
                                     a[1]*(db[0]*dcP(1,0) + db[1]*dcP(1,1) + db[2]*dcP(1,2) + db[3]*dcP(1,3))+
                                     a[2]*(db[0]*dcP(2,0) + db[1]*dcP(2,1) + db[2]*dcP(2,2) + db[3]*dcP(2,3))+
                                     a[3]*(db[0]*dcP(3,0) + db[1]*dcP(3,1) + db[2]*dcP(3,2) + db[3]*dcP(3,3)));
  secDerivs(2,2) = dzInv * dzInv *
                   (a[0]*(b[0]*dot(d2c, p00)+b[1]*dot(d2c, p01) + b[2]*dot(d2c, p02) + b[3]*dot(d2c, p03))+
                    a[1]*(b[0]*dot(d2c, p10)+b[1]*dot(d2c, p11) + b[2]*dot(d2c, p12) + b[3]*dot(d2c, p13))+
                    a[2]*(b[0]*dot(d2c, p20)+b[1]*dot(d2c, p21) + b[2]*dot(d2c, p22) + b[3]*dot(d2c, p23))+
                    a[3]*(b[0]*dot(d2c, p30)+b[1]*dot(d2c, p31) + b[2]*dot(d2c, p32) + b[3]*dot(d2c, p33)));
  return dot(a,bcP);
}

#include "Numerics/BsplineOneDimSolvers.h"

template<typename T, bool PBC0, bool PBC1, bool PBC2>
void
TricubicBsplineGridBC<T,PBC0,PBC1,PBC2>::SolvePeriodicInterp(const Array<T,3>& data, Array<T,3>& P)
{
  std::vector<T> dTemp(Nx),pTemp(Nx);
  // Do X direction
  for (int iy=0; iy<Ny; iy++)
    for (int iz=0; iz<Nz; iz++)
    {
      for(int ix=0; ix<Nx; ix++)
        dTemp[ix]=data(ix,iy,iz);
      SolvePeriodicInterp1D<T>::apply(dTemp,pTemp);
      for(int ix=0; ix<Nx; ix++)
        P(ix+1,iy+1,iz+1)=pTemp[ix];
      //SolvePeriodicInterp1D(data(Range(0,Nx-1), iy, iz), P(Range(1,Nx),iy+1, iz+1));
    }
  dTemp.resize(Ny);
  pTemp.resize(Ny);
  // Do Y direction
  for (int ix=0; ix<Nx; ix++)
    for (int iz=0; iz<Nz; iz++)
    {
      for(int iy=0; iy<Ny; iy++)
        dTemp[iy]=P(ix+1,iy+1,iz+1);
      SolvePeriodicInterp1D<T>::apply(dTemp,pTemp);
      for(int iy=0; iy<Ny; iy++)
        P(ix+1,iy+1,iz+1)=pTemp[iy];
      //SolvePeriodicInterp1D(P(ix+1,Range(1,Ny), iz+1), P(ix+1, Range(1,Ny), iz+1));
    }
  dTemp.resize(Nz);
  pTemp.resize(Nz);
  // Do z direction
  for (int ix=0; ix<Nx; ix++)
    for (int iy=0; iy<Ny; iy++)
    {
      for(int iz=0; iz<Nz; iz++)
        dTemp[iz]=P(ix+1,iy+1,iz+1);
      SolvePeriodicInterp1D<T>::apply(dTemp,pTemp);
      for(int iz=0; iz<Nz; iz++)
        P(ix+1,iy+1,iz+1)=pTemp[iz];
      //SolvePeriodicInterp1D(P(ix+1,iy+1,Range(1,Nz)), P(ix+1, iy+1, Range(1,Nz)));
    }
}

template<typename T, bool PBC0, bool PBC1, bool PBC2>
void
TricubicBsplineGridBC<T,PBC0,PBC1,PBC2>::MakePeriodic(Array<T,3>& P)
{
  // Now, make periodic
  for (int ix=0; ix<(Nx+3); ix++)
    for (int iy=0; iy<(Ny+3); iy++)
      for (int iz=0; iz<(Nz+3); iz++)
        P(ix, iy, iz) = P((ix+Nx-1)%Nx+1, (iy+Ny-1)%Ny+1, (iz+Nz-1)%Nz+1);
}

template<typename T, bool PBC0, bool PBC1, bool PBC2>
void
TricubicBsplineGridBC<T,PBC0,PBC1,PBC2>::SolveFirstDerivInterp(const Array<T,3>& data, Array<T,3>& P)
{
  real_type bcLower[]= {-3.0,0.0,3.0,0.0};
  real_type bcUpper[]= {-3.0,0.0,3.0,0.0};
  int Mx=Nx+2;
  int My=Ny+2;
  int Mz=Nz+2;
  // Do X direction
  std::vector<T> dTemp(Mx),pTemp(Mx);
  for (int iy=0; iy<Ny; iy++)
    for (int iz=0; iz<Nz; iz++)
    {
      for(int ix=0; ix<Nx; ix++)
        dTemp[ix]=data(ix,iy,iz);
      SolveFirstDerivInterp1D<T>::apply(dTemp,pTemp,Nx,bcLower,bcUpper);
      for(int ix=0; ix<Nx; ix++)
        P(ix,iy+1,iz+1)=pTemp[ix];
      //SolveDerivInterp1D(data(Range::all(), iy, iz),  P(Range::all(), iy+1, iz+1), BC, BC);
    }
  // Do Y direction
  dTemp.resize(My);
  pTemp.resize(My);
  for (int ix=0; ix<Mx; ix++)
    for (int iz=0; iz<Mz; iz++)
    {
      for(int iy=0; iy<Ny; iy++)
        dTemp[iy]=P(ix,iy+1,iz);
      SolveFirstDerivInterp1D<T>::apply(dTemp,pTemp,Ny,bcLower,bcUpper);
      for(int iy=0; iy<My; iy++)
        P(ix,iy,iz)=pTemp[iy];
      //SolveDerivInterp1D(P(ix,Range(1,Ny), iz),  P(ix, Range::all(), iz), BC, BC);
    }
  // Do z direction
  dTemp.resize(Mz);
  pTemp.resize(Mz);
  for (int ix=0; ix<Mx; ix++)
    for (int iy=0; iy<My; iy++)
    {
      for(int iz=0; iz<Nz; iz++)
        dTemp[iz]=P(ix,iy,iz+1);
      SolveFirstDerivInterp1D<T>::apply(dTemp,pTemp,Nz,bcLower,bcUpper);
      for(int iz=0; iz<Mz; iz++)
        P(ix,iy,iz)=pTemp[iz];
      //SolveDerivInterp1D(P(ix, iy, Range(1,Nz)),  P(ix, iy, Range::all()), BC, BC);
    }
}

template<typename T, bool PBC0, bool PBC1, bool PBC2>
void
TricubicBsplineGridBC<T,PBC0,PBC1,PBC2>::Init(const Array<T,3>& data, Array<T,3>& P)
{
  //if (Periodic) {
  if(PBC0||PBC1||PBC2) // if any direction is periodic, make it periodic
  {
    P.resize(Nx+3,Ny+3,Nz+3);
    if (Interpolating)
      SolvePeriodicInterp(data,P);
    else
    {
      for (int ix=0; ix<Nx; ix++)
        for (int iy=0; iy<Ny; iy++)
          for (int iz=0; iz<Nz; iz++)
            P(ix+1,iy+1,iz+1)=data(ix,iy,iz);
      //P(Range(1,Nx),Range(1,Ny),Range(1,Nz)) = data;
    }
    MakePeriodic(P);
  }
  else
  {
    P.resize(Nx+2,Ny+2,Nz+2);
    SolveFirstDerivInterp(data,P);
    //cerr << "Nonperiodic Tricubic B-splines not yet supported.\n";
    //abort();
  }
}
