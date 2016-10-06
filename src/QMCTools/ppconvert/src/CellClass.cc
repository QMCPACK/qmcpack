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
    
    



#include "CellClass.h"
#include "common/Communication.h"

void
CellClass::Broadcast (CommunicatorClass &comm, int root)
{
  Mat3 lattice = Lattice.GetDirect();
  comm.Broadcast (root, lattice);
  Lattice.SetDirect(lattice);
  int numIons = IonPos.size();
  comm.Broadcast (root, numIons);
  IonPos.resize(numIons);
  AtomTypes.resize(numIons);
  comm.Broadcast (root, IonPos);
  comm.Broadcast (root, AtomTypes);

  // Now send the GVecs
  GVecs.Broadcast (comm, root);

  // Setup the FFT
  FFT.Setup();
}


void
CellClass::SetLattice(Mat3 A)
{
  Lattice.SetDirect(A);
}
