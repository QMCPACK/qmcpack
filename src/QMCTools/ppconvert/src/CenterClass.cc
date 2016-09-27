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
    
    



#include "CenterClass.h"
#include "CellClass.h"
#include "common/Communication.h"

void
CenterClass::Broadcast(CommunicatorClass &comm, int root)
{
  comm.Broadcast (root, r);
  comm.Broadcast (root, Radius);
  comm.Broadcast (root, Spherical);
  comm.Broadcast (root, NumOrbitals);
  comm.Broadcast (root, NumUp);
  comm.Broadcast (root, NumDown);
  int numSites = IdenticalSites.size();
  comm.Broadcast (root, numSites);
  if (comm.MyProc() != root)
    IdenticalSites.resize(numSites);
  comm.Broadcast (root, IdenticalSites);
}


void
CenterClass::SetupBitfield (CellClass &cell)
{
  assert (sizeof(int) == 4);
  FFTBox &FFT = cell.FFT;
  cell.SetupFFT();

  int nx, ny, nz;
  FFT.GetDims(nx,ny,nz);
  double sx, sy, sz;
  double nxInv = 1.0/(double)nx;
  double nyInv = 1.0/(double)ny;
  double nzInv = 1.0/(double)nz;
  Vec3 r2;
  Vec3 a0 = cell.Lattice.a(0);
  Vec3 a1 = cell.Lattice.a(1);
  Vec3 a2 = cell.Lattice.a(2);
  Bitfield.Init (nx, ny, nz);
  Indices.clear();
  Offsets.clear();
  for (int ix=0; ix<nx; ix++) {
    sx = (double)ix * nxInv;
    for (int iy=0; iy<ny; iy++) {
      sy = (double)iy * nyInv;
      for (int iz=0; iz<nz; iz++) {
	sz = (double) iz * nzInv;
	r2 = sx*a0 + sy*a1 + sz*a2;
	Vec3 diff = cell.Lattice.MinImage(r - r2);
	Bitfield.Set (ix,iy,iz, (dot(diff,diff) < Radius*Radius));
	if (dot (diff,diff) < Radius*Radius) {
	  Indices.push_back(Int3(ix,iy,iz));
	  int offset = &(FFT.rBox(ix,iy,iz)) - &(FFT.rBox(0,0,0));
	  Offsets.push_back(offset);
	}
      }
    }
  }
}
