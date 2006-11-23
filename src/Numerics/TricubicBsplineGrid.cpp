//////////////////////////////////////////////////////////////////
// (c) Copyright 2006-  Kenneth Esler 
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

template<typename T>
inline TricubicBsplineGrid<T>::TricubicBsplineGrid()
{
  A(0,0) = -1.0/6.0; A(0,1) =  3.0/6.0; A(0,2) = -3.0/6.0; A(0,3) = 1.0/6.0;
  A(1,0) =  3.0/6.0; A(1,1) = -6.0/6.0; A(1,2) =  3.0/6.0; A(1,3) = 0.0/6.0;
  A(2,0) = -3.0/6.0; A(2,1) =  0.0/6.0; A(2,2) =  3.0/6.0; A(2,3) = 0.0/6.0;
  A(3,0) =  1.0/6.0; A(3,1) =  4.0/6.0; A(3,2) =  1.0/6.0; A(3,3) = 0.0/6.0;

  dA(0,0)= 0.0; dA(0,1)= 0.0; dA(0,2)= 0.0; dA(0,3)= 0.0;
  dA(1,0)=-0.5; dA(1,1)= 1.5; dA(1,2)=-1.5; dA(1,3)= 0.5;
  dA(2,0)= 1.0; dA(2,1)=-2.0; dA(2,2)= 1.0; dA(2,3)= 0.0;
  dA(3,0)=-0.5; dA(3,1)= 0.0; dA(3,2)= 0.5; dA(3,3)= 0.0;

  d2A(0,0)= 0.0; d2A(0,1)= 0.0; d2A(0,2)= 0.0; d2A(0,3)= 0.0;
  d2A(1,0)= 0.0; d2A(1,1)= 0.0; d2A(1,2)= 0.0; d2A(1,3)= 0.0;
  d2A(2,0)=-1.0; d2A(2,1)= 3.0; d2A(2,2)=-3.0; d2A(2,3)= 1.0;
  d2A(3,0)= 1.0; d2A(3,1)=-2.0; d2A(3,2)= 1.0; d2A(3,3)= 0.0;

  d3A(0,0)= 0.0; d3A(0,1)= 0.0; d3A(0,2)= 0.0; d3A(0,3)= 0.0;
  d3A(1,0)= 0.0; d3A(1,1)= 0.0; d3A(1,2)= 0.0; d3A(1,3)= 0.0;
  d3A(2,0)= 0.0; d3A(2,1)= 0.0; d3A(1,2)= 2.0; d3A(2,3)= 0.0;
  d3A(3,0)=-1.0; d3A(3,1)= 3.0; d3A(3,2)=-3.0; d3A(3,3)= 1.0;
}

template<typename T>
inline TricubicBsplineGrid<T>::TricubicBsplineGrid(const TricubicBsplineGrid<T>& rhs)
{
  setGrid(rhs.xStart,rhs.xEnd,rhs.yStart,rhs.yEnd,rhs.zStart,rhs.zEnd,
      rhs.Nx, rhs.Ny, rhs.Nz, rhs.Interpolating, rhs.Periodic);
}

template<typename T>
inline TricubicBsplineGrid<T>&
TricubicBsplineGrid<T>::operator=(const TricubicBsplineGrid<T>& rhs)
{
  setGrid(rhs.xStart,rhs.xEnd,rhs.yStart,rhs.yEnd,rhs.zStart,rhs.zEnd,
      rhs.Nx, rhs.Ny, rhs.Nz, rhs.Interpolating, rhs.Periodic);
  return *this;
}

template<typename T>
inline void TricubicBsplineGrid<T>::setGrid(real_type xi, real_type xf, 
    real_type yi, real_type yf, real_type zi, real_type zf, 
    int nx, int ny, int nz, bool interp, bool periodic, bool openend)
{
  Interpolating=interp;
  Periodic=periodic;
  xStart=xi;  xEnd=xf;  yStart=yi;  yEnd=yf;  zStart=zi; zEnd=zf;
  Lx    = xf-xi;  Ly    = yf-yi;  Lz    = zf-zi;
  LxInv = 1.0/Lx; LyInv = 1.0/Ly; LzInv = 1.0/Lz;
  if(openend) {
    Nx=nx; Ny=ny; Nz=nz;
  } else {
    Nx=nx-1; Ny=ny-1; Nz=nz-1;
  }
  dx = Lx/static_cast<real_type>(Nx); dxInv = 1.0/dx;
  dy = Ly/static_cast<real_type>(Ny); dyInv = 1.0/dy;
  dz = Lz/static_cast<real_type>(Nz); dzInv = 1.0/dz;
}

template<typename T>
inline void TricubicBsplineGrid<T>::Find(real_type x, real_type y, real_type z)
{
  real_type xDelta = x - xStart;
  real_type yDelta = y - yStart;
  real_type zDelta = z - zStart;

  if (Periodic) {
    //     xDelta -= nearbyint(xDelta*LxInv)*Lx;
    //     yDelta -= nearbyint(yDelta*LyInv)*Ly;
    //     zDelta -= nearbyint(zDelta*LzInv)*Lz;
    xDelta -= std::floor(xDelta*LxInv)*Lx;
    yDelta -= std::floor(yDelta*LyInv)*Ly;
    zDelta -= std::floor(zDelta*LzInv)*Lz;
  }

  real_type xInt, yInt, zInt;
  real_type tx = modf (xDelta*dxInv, &xInt);
  real_type ty = modf (yDelta*dyInv, &yInt);
  real_type tz = modf (zDelta*dzInv, &zInt);
  int ix = (int)xInt;
  int iy = (int)yInt;
  int iz = (int)zInt;
  ix0=ix; ix1=ix+1; ix2=ix+2; ix3=ix+3;
  iy0=iy; iy1=iy+1; iy2=iy+2; iy3=iy+3;
  iz0=iz; iz1=iz+1; iz2=iz+2; iz3=iz+3;

  px[0] = tx*tx*tx;  py[0] = ty*ty*ty; pz[0] = tz*tz*tz;
  px[1] = tx*tx;     py[1] = ty*ty;    pz[1] = tz*tz;
  px[2] = tx;        py[2] = ty;       pz[2] = tz;
  px[3] = 1.0;       py[3] = 1.0;      pz[3] = 1.0;

  a=dot(px,A);
  b=dot(py,A);
  c=dot(pz,A);
}

template<typename T>
inline void TricubicBsplineGrid<T>::FindAll(real_type x, real_type y, real_type z) {
  Find(x,y,z);
  //da=dot(px,dA);
  //db=dot(py,dA);
  //dc=dot(pz,dA);
  //// First derivatives;
  da[0] = dA(1,0)*px[1]; da[1] =dA(1,1)*px[1]; da[2] =dA(1,2)*px[1]; da[3] =dA(1,3)*px[1];
  da[0]+= dA(2,0)*px[2]; da[1]+=dA(2,1)*px[2]; da[2]+=dA(2,2)*px[2]; da[3]+=dA(2,3)*px[2];
  da[0]+= dA(3,0)*px[3]; da[1]+=dA(3,1)*px[3]; da[2]+=dA(3,2)*px[3]; da[3]+=dA(3,3)*px[3];  

  db[0] = dA(1,0)*py[1]; db[1] =dA(1,1)*py[1]; db[2] =dA(1,2)*py[1]; db[3] =dA(1,3)*py[1];
  db[0]+= dA(2,0)*py[2]; db[1]+=dA(2,1)*py[2]; db[2]+=dA(2,2)*py[2]; db[3]+=dA(2,3)*py[2];
  db[0]+= dA(3,0)*py[3]; db[1]+=dA(3,1)*py[3]; db[2]+=dA(3,2)*py[3]; db[3]+=dA(3,3)*py[3];  

  dc[0] = dA(1,0)*pz[1]; dc[1] =dA(1,1)*pz[1]; dc[2] =dA(1,2)*pz[1]; dc[3] =dA(1,3)*pz[1];
  dc[0]+= dA(2,0)*pz[2]; dc[1]+=dA(2,1)*pz[2]; dc[2]+=dA(2,2)*pz[2]; dc[3]+=dA(2,3)*pz[2];
  dc[0]+= dA(3,0)*pz[3]; dc[1]+=dA(3,1)*pz[3]; dc[2]+=dA(3,2)*pz[3]; dc[3]+=dA(3,3)*pz[3];  

  // Second derivatives
  //d2a=dot(px,d2A);
  //d2b=dot(py,d2A);
  //d2c=dot(pz,d2A);
  d2a[0] = d2A(2,0)*px[2]; d2a[1] =d2A(2,1)*px[2]; d2a[2] =d2A(2,2)*px[2]; d2a[3] =d2A(2,3)*px[2];
  d2a[0]+= d2A(3,0)*px[3]; d2a[1]+=d2A(3,1)*px[3]; d2a[2]+=d2A(3,2)*px[3]; d2a[3]+=d2A(3,3)*px[3];  

  d2b[0] = d2A(2,0)*py[2]; d2b[1] =d2A(2,1)*py[2]; d2b[2] =d2A(2,2)*py[2]; d2b[3] =d2A(2,3)*py[2];
  d2b[0]+= d2A(3,0)*py[3]; d2b[1]+=d2A(3,1)*py[3]; d2b[2]+=d2A(3,2)*py[3]; d2b[3]+=d2A(3,3)*py[3];  

  d2c[0] = d2A(2,0)*pz[2]; d2c[1] =d2A(2,1)*pz[2]; d2c[2] =d2A(2,2)*pz[2]; d2c[3] =d2A(2,3)*pz[2];
  d2c[0]+= d2A(3,0)*pz[3]; d2c[1]+=d2A(3,1)*pz[3]; d2c[2]+=d2A(3,2)*pz[3]; d2c[3]+=d2A(3,3)*pz[3];  
}

template<typename T>
inline T TricubicBsplineGrid<T>::evaluate(const Array<T,3>& P) const
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

template<typename T>
inline T TricubicBsplineGrid<T>::evaluate(const Array<T,3>& P,
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

template<typename T> 
inline T TricubicBsplineGrid<T>::evaluate(const Array<T,3>& P,
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
  TinyMatrix<T,4,4> cP, dcP;
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

template<typename T> void
TricubicBsplineGrid<T>::SolvePeriodicInterp(const Array<T,3>& data, Array<T,3>& P)
{
  vector<T> dTemp(Nx),pTemp(Nx);
  // Do X direction
  for (int iy=0; iy<Ny; iy++)
    for (int iz=0; iz<Nz; iz++)  
    {
      for(int ix=0; ix<Nx;ix++) dTemp[ix]=data(ix,iy,iz);
      SolvePeriodicInterp1D<T>::apply(dTemp,pTemp);
      for(int ix=0; ix<Nx;ix++) P(ix+1,iy+1,iz+1)=pTemp[ix];
      //SolvePeriodicInterp1D(data(Range(0,Nx-1), iy, iz), P(Range(1,Nx),iy+1, iz+1));
    }
  

  dTemp.resize(Ny);
  pTemp.resize(Ny);
  // Do Y direction
  for (int ix=0; ix<Nx; ix++)
    for (int iz=0; iz<Nz; iz++) 
    {
      for(int iy=0;iy<Ny; iy++) dTemp[iy]=P(ix+1,iy+1,iz+1);
      SolvePeriodicInterp1D<T>::apply(dTemp,pTemp);
      for(int iy=0;iy<Ny; iy++) P(ix+1,iy+1,iz+1)=pTemp[iy];
      //SolvePeriodicInterp1D(P(ix+1,Range(1,Ny), iz+1), P(ix+1, Range(1,Ny), iz+1));
    }
  
  dTemp.resize(Nz);
  pTemp.resize(Nz);
  // Do z direction
  for (int ix=0; ix<Nx; ix++)
    for (int iy=0; iy<Ny; iy++) 
    {
      for(int iz=0; iz<Nz; iz++) dTemp[iz]=P(ix+1,iy+1,iz+1);
      SolvePeriodicInterp1D<T>::apply(dTemp,pTemp);
      for(int iz=0; iz<Nz; iz++) P(ix+1,iy+1,iz+1)=pTemp[iz];
      //SolvePeriodicInterp1D(P(ix+1,iy+1,Range(1,Nz)), P(ix+1, iy+1, Range(1,Nz)));
    }

}

template<typename T> void
TricubicBsplineGrid<T>::MakePeriodic(Array<T,3>& P)
{
  // Now, make periodic
  for (int ix=0; ix<(Nx+3); ix++)
    for (int iy=0; iy<(Ny+3); iy++)
      for (int iz=0; iz<(Nz+3); iz++)
	P(ix, iy, iz) = P((ix+Nx-1)%Nx+1, (iy+Ny-1)%Ny+1, (iz+Nz-1)%Nz+1);
}

template<typename T> void
TricubicBsplineGrid<T>::Init(const Array<T,3>& data, Array<T,3>& P)
{
  P.resize(Nx+3,Ny+3,Nz+3);
  if (Periodic) {
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
  else {
    cerr << "Nonperiodic Tricubic B-splines not yet supported.\n";
    abort();
  }
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/
