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
    
    



#ifndef CENTER_CLASS_H
#define CENTER_CLASS_H

#include "common/Blitz.h"
#include <vector>

typedef TinyVector<int,3> Int3;

class CellClass;
class CommunicatorClass;

class BitField3DClass
{
  int Nx, Ny, Nz;
  Array<unsigned int,1> Field;
public:

  inline bool operator()(int ix, int iy, int iz) {
    int num = iz + Nz*(iy + Ny*ix);
    int wordnum = num >> 5;
    int bitnum  = num & 0x001f;
    unsigned int mask = 1 << bitnum;
    return (Field(wordnum) & mask) ? true : false;
  }

  inline double GetDouble(int ix, int iy, int iz) {
    int num = iz + Nz*(iy + Ny*ix);
    int wordnum = num >> 5;
    int bitnum  = num & 0x001f;
    unsigned int mask = 1 << bitnum;
    return (Field(wordnum) & mask) ? 1.0 : 0.0;
  }
  
  inline void Set (int ix, int iy, int iz, bool val)
  {
    int num = iz + Nz*(iy + Ny*ix);
    int wordnum = num >> 5;
    int bitnum  = num & 0x001f;
    unsigned int mask = 1 << bitnum;
    if (val)
      Field(wordnum) |= mask;
    else 
      Field(wordnum) &= (mask ^ 0xffffffff);
  }
  inline void Init (int nx, int ny, int nz) 
  {
    Nx = nx; Ny = ny; Nz = nz;
    Field.resize((nx*ny*nz+31)/32);
    Field = 0;
  }
};

class CenterClass
{
public:
  Vec3 r;
  // This stores a list of sites which are identical, by symmetry, to
  // the one contained in r.  This array contains r for convenience.
  Array<Vec3,1> IdenticalSites;
  // This stores whether to reflect each orbital along each of the 3
  // planes, inverting the size of the X, Y, or Z coordinate.  Each
  // value should be 1.0 or -1.0;
  Array<Vec3,1> Reflections;
  double Radius, SkinThickness;
  bool Spherical;
  int NumOrbitals;
  int NumUp, NumDown;
  BitField3DClass Bitfield;
  void Broadcast (CommunicatorClass &comm, int root=0);
  std::vector<Int3> Indices;
  std::vector<int> Offsets;

  void SetupBitfield (CellClass &cell);

  CenterClass (Vec3 center, double radius, 
	       int numOrb, bool spherical=true)
  {
    r = center;
    Radius = radius;
    Spherical = spherical;
    NumOrbitals = numOrb;
  }
  CenterClass() : 
    Spherical(true), r(0.0, 0.0, 0.0), Radius(0.0), NumOrbitals(1),
    SkinThickness(0.0)
  {

  }
};

#endif
