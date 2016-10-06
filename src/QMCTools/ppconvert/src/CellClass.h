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
    
    



#ifndef CELL_CLASS_H
#define CELL_CLASS_H

#include "LatticeClass.h"
#include "common/FFTBox.h"

class CellClass
{
private:
  bool FFTisSetup;
public:
  LatticeClass Lattice;
  GVecsClass GVecs;
  FFTBox FFT;
  Array<Vec3, 1> IonPos;
  Array<Vec3, 1> IonForces;
  Array<int,  1> AtomTypes;
  Array<int,  1> Zion;

  /////////////////////
  // Member functions//
  /////////////////////
  void SetLattice (Mat3 A);
  inline void SetupFFT()
  {  
    if (!FFTisSetup) {
      FFTisSetup = true;
      FFT.Setup();
    }
  }  

  void Broadcast (CommunicatorClass &comm, int root);

  // Constructor
  CellClass() : FFT(GVecs), FFTisSetup(false)
  {
  }
};

#endif
