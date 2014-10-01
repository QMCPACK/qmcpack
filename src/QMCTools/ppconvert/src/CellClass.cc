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
