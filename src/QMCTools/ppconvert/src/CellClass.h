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
