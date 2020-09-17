
#ifndef QMCPLUSPLUS_AFQMC_ONEBODYHAMILTONIAN_H
#define QMCPLUSPLUS_AFQMC_ONEBODYHAMILTONIAN_H

#include <iostream>
#include <vector>
#include <map>
#include <fstream>

#include "AFQMC/config.h"

namespace qmcplusplus
{
namespace afqmc
{
class OneBodyHamiltonian : public AFQMCInfo
{
public:
  OneBodyHamiltonian(const AFQMCInfo& info, boost::multi::array<ValueType, 2>&& h, ValueType nuc = 0, ValueType frz = 0)
      : AFQMCInfo(info), H1(h), NuclearCoulombEnergy(nuc), FrozenCoreEnergy(frz)
  {
    name = "OneBodyHamiltonian";
  }

  OneBodyHamiltonian(OneBodyHamiltonian const& other) = delete;
  OneBodyHamiltonian(OneBodyHamiltonian&& other)      = default;
  OneBodyHamiltonian& operator=(OneBodyHamiltonian const& other) = delete;
  OneBodyHamiltonian& operator=(OneBodyHamiltonian&& other) = delete;

  ~OneBodyHamiltonian() {}

  boost::multi::array<ValueType, 2> getH1() const
  {
    boost::multi::array<ValueType, 2> H_(H1);
    return H_;
  }

  // this should never be used outside initialization routines.
  ValueType H(IndexType I, IndexType J) const
  {
    if ((I >= NMO && J < NMO) || (I < NMO && J >= NMO))
      return ValueType(0);
    I = (I >= NMO) ? (I - NMO) : (I);
    J = (J >= NMO) ? (J - NMO) : (J);
    //return ValueType(0.0);
    return ValueType(H1[I][J]);
  }

  ValueType getNuclearCoulombEnergy() const { return NuclearCoulombEnergy; }

protected:
  // this should be in shared memory, no need for all the copies!!!
  boost::multi::array<ValueType, 2> H1;

  // nuclear coulomb term
  ValueType NuclearCoulombEnergy;
  ValueType FrozenCoreEnergy;
};
} // namespace afqmc
} // namespace qmcplusplus
#endif
