/////////////////////////////////////////////////////////////
// Copyright (C) 2003-2006 Bryan Clark and Kenneth Esler   //
//                                                         //
// This program is free software; you can redistribute it  //
// and/or modify it under the terms of the GNU General     //
// Public License as published by the Free Software        //
// Foundation; either version 2 of the License, or         //
// (at your option) any later version.  This program is    //
// distributed in the hope that it will be useful, but     //
// WITHOUT ANY WARRANTY; without even the implied warranty //
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. //  
// See the GNU General Public License for more details.    //
// For more information, please see the PIMC++ Home Page:  //
//           http://pathintegrals.info                     //
/////////////////////////////////////////////////////////////

#ifndef FFT_BOX_H
#define FFT_BOX_H

#include "GVecs.h"
#include "FFT.h"

class FFTBox : public FFT3D
{
protected:
  int Nx, Ny, Nz;
public:
  // The G-vectors associated with this FFT
  GVecsClass &GVecs;

  // Puts a linear c-vector into the fft box.
  void PutkVec (const zVec &c);
  
  // Puts appropriate elements of fftbox into a c vector.
  void GetkVec (zVec &c);
  
  // Adds the appropriate elements of the fftbox to a c vector.
  void AddFromVec (const zVec &c);

  // Add the c vector to the contents of the fftbox.
  void AddToVec (zVec &c);
  
  // Get FFT box dimensions
  void GetDims(int &nx, int &ny, int &nz)
  { nx=Nx; ny=Ny; nz=Nz; }

  void Setup(double factor=1.0)
  {
    GVecs.GetFFTBoxSize(Nx, Ny, Nz);
    resize(Nx, Ny, Nz);
  }

  FFTBox& operator=(const FFTBox& fft);
  FFTBox (GVecsClass &gvecs) : GVecs (gvecs)
  {
  }
};


class FFTVecBox : public FFTVec3D
{
protected:
  GVecsClass &GVecs;
  int Nx, Ny, Nz;
public:
  // Puts a linear c-vector into the fft box.
  void PutkVec (const zVecVec &c);
  
  // Puts appropriate elements of fftbox into a c vector.
  void GetkVec (zVecVec &c);
  
  // Adds the appropriate elements of the fftbox to a c vector.
  void AddFromVec (const zVecVec &c);

  // Add the c vector to the contents of the fftbox.
  void AddToVec (zVecVec &c);
  
  // Get FFT box dimensions
  void GetDims(int &nx, int &ny, int &nz)
  { nx=Nx; ny=Ny; nz=Nz; }

  void Setup()
  {
    GVecs.GetFFTBoxSize(Nx, Ny, Nz);
    resize(Nx, Ny, Nz);
  }

  FFTVecBox (GVecsClass &gvecs) : GVecs (gvecs)
  {
  }
};

inline
TinyMatrix<complex<float>,3,3> conv(TinyMatrix<complex<double>,3,3> val)
{
  TinyMatrix<complex<float>,3,3> v;
  v(0,0)=val(0,0); v(0,1)=val(0,1); v(0,2)=v(0,2);
  v(1,0)=val(1,0); v(1,1)=val(1,1); v(1,2)=v(1,2);
  v(2,0)=val(2,0); v(2,1)=val(2,1); v(2,2)=v(2,2);
  return (v);
}

inline
TinyMatrix<complex<double>,3,3> conv(TinyMatrix<complex<float>,3,3> val)
{
  TinyMatrix<complex<double>,3,3> v;
  v(0,0)=val(0,0); v(0,1)=val(0,1); v(0,2)=v(0,2);
  v(1,0)=val(1,0); v(1,1)=val(1,1); v(1,2)=v(1,2);
  v(2,0)=val(2,0); v(2,1)=val(2,1); v(2,2)=v(2,2);
  return (v);
}

class FFTMatBox : public FFTMat3D
{
protected:
  GVecsClass &GVecs;
  int Nx, Ny, Nz;
public:
  // Puts a linear c-vector into the fft box.
  void PutkVec (const zMatVec &c);
  
  // Puts appropriate elements of fftbox into a c vector.
  void GetkVec (zMatVec &c);
  
  // Adds the appropriate elements of the fftbox to a c vector.
  void AddFromVec (const zMatVec &c);

  // Add the c vector to the contents of the fftbox.
  void AddToVec (zMatVec &c);
  
  // Get FFT box dimensions
  void GetDims(int &nx, int &ny, int &nz)
  { nx=Nx; ny=Ny; nz=Nz; }

  void Setup()
  {
    GVecs.GetFFTBoxSize(Nx, Ny, Nz);
    resize(Nx, Ny, Nz);
  }

  FFTMatBox (GVecsClass &gvecs) : GVecs (gvecs)
  {
  }
};
#endif
