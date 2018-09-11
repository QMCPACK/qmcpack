
#include "QMCWaveFunctions/Fermion/DiracDeterminantSingle.h"

namespace qmcplusplus
{
DiracDeterminantSingle::DiracDeterminantSingle(SPOSetPtr const &spos, int first):
  DiracDeterminantBase(first), Phi(spos)
{
  Optimizable=false;
  if(Phi->Optimizable)
    Optimizable=true;
  OrbitalName="DiracDeterminantSingle";
  registerTimers();
}

}
