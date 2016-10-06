//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Paul R. C. Kent, kentpr@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Paul R. C. Kent, kentpr@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
    
    
//           http://pathintegrals.info                     //
/////////////////////////////////////////////////////////////

#include "FFTBox.h"

void
FFTBox::PutkVec (const zVec &c)
{
  kBox = 0.0;
  if (c.size() == GVecs.size()) 
    for (int n=0; n<c.size(); n++) 
      kBox(GVecs.Index(n)) = c(n);
  else if (c.size() == GVecs.DeltaSize()) 
    for (int n=0; n<c.size(); n++) 
      kBox(GVecs.DeltaI(n)) = c(n);
  else {
    cerr << "Incommensurate dimensions in PutkVec.\n";
    abort();
  }
}


void 
FFTBox::GetkVec (zVec &c)
{
  if (c.size() == GVecs.size()) 
    for (int n=0; n<c.size(); n++) 
      c(n) = kBox(GVecs.Index(n));
  else if (c.size() == GVecs.DeltaSize()) 
    for (int n=0; n<c.size(); n++) 
      c(n) = kBox(GVecs.DeltaI(n));
  else {
    cerr << "Incommensurate dimensions in GetkVec.\n";
    abort();
  }
}

void
FFTBox::AddFromVec (const zVec &c)
{
  if (c.size() == GVecs.size()) 
    for (int n=0; n<c.size(); n++) 
      kBox(GVecs.Index(n)) += c(n);
  else if (c.size() == GVecs.DeltaSize()) 
    for (int n=0; n<c.size(); n++) 
      kBox(GVecs.DeltaI(n)) += c(n);
  else {
    cerr << "Incommensurate dimensions in AddFromVec.\n";
    abort();
  }
}

void 
FFTBox::AddToVec (zVec &c)
{
  if (c.size() == GVecs.size()) 
    for (int n=0; n<c.size(); n++) 
      c(n) += kBox(GVecs.Index(n));
  else if (c.size() == GVecs.DeltaSize()) 
    for (int n=0; n<c.size(); n++) 
      c(n) += kBox(GVecs.DeltaI(n));
  else {
    cerr << "Incommensurate dimensions in GetkVec.\n";
    abort();
  }
}


////////////////////////////////
// FFTVecBox Member Functions //
////////////////////////////////
void
FFTVecBox::PutkVec (const zVecVec &c)
{
  complex<double> zero(0.0, 0.0);
  cVec3 zero3(zero, zero, zero);
  kBox = zero3;
  if (c.size() == GVecs.size()) 
    for (int n=0; n<c.size(); n++) {
//       int i = (GVecs.Index(n)[0]+Nx)%Nx;
//       int j = (GVecs.Index(n)[1]+Ny)%Ny;
//       int k = (GVecs.Index(n)[2]+Nz)%Nz;
//       kBox(i,j,k) = c(n);
      kBox(GVecs.Index(n)) = c(n);
    }
  else if (c.size() == GVecs.DeltaSize()) 
    for (int n=0; n<c.size(); n++) {
//       int i = (GVecs.DeltaI(n)[0]+Nx)%Nx;
//       int j = (GVecs.DeltaI(n)[1]+Ny)%Ny;
//       int k = (GVecs.DeltaI(n)[2]+Nz)%Nz;
//       kBox(i,j,k) = c(n);
      kBox(GVecs.DeltaI(n)) = c(n);
    }
  else {
    cerr << "Incommensurate dimensions in PutkVec.\n";
    abort();
  }
}

void 
FFTVecBox::GetkVec (zVecVec &c)
{
  if (c.size() == GVecs.size()) 
    for (int n=0; n<c.size(); n++) {
//       int i = (GVecs.Index(n)[0]+Nx)%Nx;
//       int j = (GVecs.Index(n)[1]+Ny)%Ny;
//       int k = (GVecs.Index(n)[2]+Nz)%Nz;
//       c(n) = kBox(i,j,k);
      c(n) = kBox(GVecs.Index(n));
    }
  else if (c.size() == GVecs.DeltaSize()) 
    for (int n=0; n<c.size(); n++) {
//       int i = (GVecs.DeltaI(n)[0]+Nx)%Nx;
//       int j = (GVecs.DeltaI(n)[1]+Ny)%Ny;
//       int k = (GVecs.DeltaI(n)[2]+Nz)%Nz;
//       c(n) = kBox(i,j,k);
      c(n) = kBox(GVecs.DeltaI(n));
    }
  else {
    cerr << "Incommensurate dimensions in GetkVec.\n";
    abort();
  }
}


void
FFTVecBox::AddFromVec (const zVecVec &c)
{
  if (c.size() == GVecs.size()) 
    for (int n=0; n<c.size(); n++) {
//       int i = (GVecs.Index(n)[0]+Nx)%Nx;
//       int j = (GVecs.Index(n)[1]+Ny)%Ny;
//       int k = (GVecs.Index(n)[2]+Nz)%Nz;
//       kBox(i,j,k) += c(n);
      kBox(GVecs.Index(n)) += c(n);
    }
  else if (c.size() == GVecs.DeltaSize()) 
    for (int n=0; n<c.size(); n++) {
//       int i = (GVecs.DeltaI(n)[0]+Nx)%Nx;
//       int j = (GVecs.DeltaI(n)[1]+Ny)%Ny;
//       int k = (GVecs.DeltaI(n)[2]+Nz)%Nz;
//       kBox(i,j,k) += c(n);
      kBox(GVecs.DeltaI(n)) += c(n);
    }
  else {
    cerr << "Incommensurate dimensions in AddFromVec.\n";
    abort();
  }
}

void 
FFTVecBox::AddToVec (zVecVec &c)
{
  if (c.size() == GVecs.size()) 
    for (int n=0; n<c.size(); n++) 
      c(n) += kBox(GVecs.Index(n));
  else if (c.size() == GVecs.DeltaSize()) 
    for (int n=0; n<c.size(); n++) 
      c(n) += kBox(GVecs.DeltaI(n));
  else {
    cerr << "Incommensurate dimensions in GetkVec.\n";
    abort();
  }
}


////////////////////////////////
// FFTMatBox Member Functions //
////////////////////////////////
void
FFTMatBox::PutkVec (const zMatVec &c)
{
  TinyMatrix<complex<FFT_FLOAT>,3,3> zero;
  zero(0,0)=0.0; zero(0,1)=0.0; zero(0,2)= 0.0;
  zero(1,0)=0.0; zero(1,1)=0.0; zero(1,2)= 0.0;
  zero(2,0)=0.0; zero(2,1)=0.0; zero(2,2)= 0.0;
  kBox = zero;
  if (c.size() == GVecs.size()) 
    for (int n=0; n<c.size(); n++) {
      TinyMatrix<complex<FFT_FLOAT>,3,3> &val = kBox(GVecs.Index(n));
      val(0,0)=c(n)(0,0); val(0,1)=c(n)(0,1); val(0,2)=c(n)(0,2);
      val(1,0)=c(n)(1,0); val(1,1)=c(n)(0,1); val(1,2)=c(n)(1,2);
      val(2,0)=c(n)(2,0); val(2,1)=c(n)(0,1); val(2,2)=c(n)(2,2);
    }
  else if (c.size() == GVecs.DeltaSize()) 
    for (int n=0; n<c.size(); n++) {
      TinyMatrix<complex<FFT_FLOAT>,3,3> &val = kBox(GVecs.DeltaI(n));
      val(0,0)=c(n)(0,0); val(0,1)=c(n)(0,1); val(0,2)=c(n)(0,2);
      val(1,0)=c(n)(1,0); val(1,1)=c(n)(0,1); val(1,2)=c(n)(1,2);
      val(2,0)=c(n)(2,0); val(2,1)=c(n)(0,1); val(2,2)=c(n)(2,2);
    }
  else {
    cerr << "Incommensurate dimensions in PutkVec.\n";
    abort();
  }
}

void 
FFTMatBox::GetkVec (zMatVec &c)
{
  if (c.size() == GVecs.size()) 
    for (int n=0; n<c.size(); n++) {
      TinyMatrix<complex<FFT_FLOAT>,3,3> &val = kBox(GVecs.Index(n)); 
      c(n)(0,0)=val(0,0); c(n)(0,1)=val(0,1); c(n)(0,2)=val(0,2);
      c(n)(1,0)=val(1,0); c(n)(1,1)=val(1,1); c(n)(1,2)=val(1,2);
      c(n)(2,0)=val(2,0); c(n)(2,1)=val(2,1); c(n)(2,2)=val(2,2);
    }
  else if (c.size() == GVecs.DeltaSize()) 
    for (int n=0; n<c.size(); n++) {
//       TinyMatrix<complex<FFT_FLOAT>,3,3> &val = kBox(GVecs.Index(n));
//       c(n) = kBox(GVecs.DeltaI(n));
      TinyMatrix<complex<FFT_FLOAT>,3,3> &val = kBox(GVecs.DeltaI(n)); 
      c(n)(0,0)=val(0,0); c(n)(0,1)=val(0,1); c(n)(0,2)=val(0,2);
      c(n)(1,0)=val(1,0); c(n)(1,1)=val(1,1); c(n)(1,2)=val(1,2);
      c(n)(2,0)=val(2,0); c(n)(2,1)=val(2,1); c(n)(2,2)=val(2,2);
    }
  else {
    cerr << "Incommensurate dimensions in GetkVec.\n";
    abort();
  }
}


void
FFTMatBox::AddFromVec (const zMatVec &c)
{
  if (c.size() == GVecs.size()) 
    for (int n=0; n<c.size(); n++) {
      //kBox(GVecs.Index(n)) += c(n);
      TinyMatrix<complex<FFT_FLOAT>,3,3> &val = kBox(GVecs.Index(n)); 
      val(0,0)+=c(n)(0,0); val(0,1)+=c(n)(0,1); val(0,2)+=c(n)(0,2);
      val(1,0)+=c(n)(1,0); val(1,1)+=c(n)(1,1); val(1,2)+=c(n)(1,2);
      val(2,0)+=c(n)(2,0); val(2,1)+=c(n)(2,1); val(2,2)+=c(n)(2,2);
    }
  else if (c.size() == GVecs.DeltaSize()) 
    for (int n=0; n<c.size(); n++) {
      TinyMatrix<complex<FFT_FLOAT>,3,3> &val = kBox(GVecs.DeltaI(n)); 
      val(0,0)+=c(n)(0,0); val(0,1)+=c(n)(0,1); val(0,2)+=c(n)(0,2);
      val(1,0)+=c(n)(1,0); val(1,1)+=c(n)(1,1); val(1,2)+=c(n)(1,2);
      val(2,0)+=c(n)(2,0); val(2,1)+=c(n)(2,1); val(2,2)+=c(n)(2,2);
    }
  else {
    cerr << "Incommensurate dimensions in AddFromVec.\n";
    abort();
  }
}

void 
FFTMatBox::AddToVec (zMatVec &c)
{
  if (c.size() == GVecs.size()) 
    for (int n=0; n<c.size(); n++) {
      TinyMatrix<complex<FFT_FLOAT>,3,3> &val = kBox(GVecs.Index(n)); 
      c(n)(0,0)+=val(0,0); c(n)(0,1)+=val(0,1); c(n)(0,2)+=val(0,2);
      c(n)(1,0)+=val(1,0); c(n)(1,1)+=val(1,1); c(n)(1,2)+=val(1,2);
      c(n)(2,0)+=val(2,0); c(n)(2,1)+=val(2,1); c(n)(2,2)+=val(2,2);
      // c(n) += kBox(GVecs.Index(n));
    }
  else if (c.size() == GVecs.DeltaSize()) 
    for (int n=0; n<c.size(); n++) {
      //      c(n) += kBox(GVecs.DeltaI(n));
      TinyMatrix<complex<FFT_FLOAT>,3,3> &val = kBox(GVecs.DeltaI(n)); 
      c(n)(0,0)+=val(0,0); c(n)(0,1)+=val(0,1); c(n)(0,2)+=val(0,2);
      c(n)(1,0)+=val(1,0); c(n)(1,1)+=val(1,1); c(n)(1,2)+=val(1,2);
      c(n)(2,0)+=val(2,0); c(n)(2,1)+=val(2,1); c(n)(2,2)+=val(2,2);
    }
  else {
    cerr << "Incommensurate dimensions in GetkVec.\n";
    abort();
  }
}


FFTBox&
FFTBox::operator=(const FFTBox &fft)
{
  FFT3D::operator=(fft);
  Nx = fft.Nx;
  Ny = fft.Ny;
  Nz = fft.Nz;

  return *this;
}
