
#ifndef QMCPLUSPLUS_AFQMC_REALDENSEHAMILTONIAN_V2_H
#define QMCPLUSPLUS_AFQMC_REALDENSEHAMILTONIAN_V2_H

#include <iostream>
#include <vector>
#include <map>
#include <fstream>

#include "hdf/hdf_archive.h"

#include "AFQMC/config.h"
#include "AFQMC/Memory/utilities.hpp"
#include "AFQMC/Utilities/taskgroup.h"
#include "AFQMC/Numerics/ma_operations.hpp"

#include "AFQMC/Hamiltonians/OneBodyHamiltonian.hpp"

#include "AFQMC/HamiltonianOperations/HamiltonianOperations.hpp"

namespace qmcplusplus
{
namespace afqmc
{
class RealDenseHamiltonian_v2 : public OneBodyHamiltonian
{
public:
  RealDenseHamiltonian_v2(AFQMCInfo const& info,
                          xmlNodePtr cur,
                          boost::multi::array<ValueType, 2>&& h,
                          TaskGroup_& tg_,
                          ValueType nucE = 0,
                          ValueType fzcE = 0)
      : OneBodyHamiltonian(info, std::move(h), nucE, fzcE),
        TG(tg_),
        fileName(""),
        batched("no"),
        ooc("no"),
        max_memory_MB(2000)
  {
    if (number_of_devices() > 0)
      batched = "yes";
    std::string str("yes");
    ParameterSet m_param;
    m_param.add(fileName, "filename");
    m_param.add(max_memory_MB, "max_memory");
    if (TG.TG_local().size() == 1)
      m_param.add(batched, "batched");
    if (TG.TG_local().size() == 1)
      m_param.add(ooc, "ooc");
    m_param.put(cur);

    if (omp_get_num_threads() > 1 && (batched != "yes" && batched != "true"))
    {
      app_log() << " WARNING!!!: Found OMP_NUM_THREADS > 1 with batched=no.\n"
                << "             This will lead to low performance. Set batched=yes. \n";
    }

    if (TG.TG_local().size() > 1 || not(batched == "yes" || batched == "true"))
    {
      app_log() << " alt_version = true only on GPU build with batched=true. " << std::endl;
      APP_ABORT("");
    }
  }

  ~RealDenseHamiltonian_v2() {}

  RealDenseHamiltonian_v2(RealDenseHamiltonian_v2 const& other) = delete;
  RealDenseHamiltonian_v2(RealDenseHamiltonian_v2&& other)      = default;
  RealDenseHamiltonian_v2& operator=(RealDenseHamiltonian_v2 const& other) = delete;
  RealDenseHamiltonian_v2& operator=(RealDenseHamiltonian_v2&& other) = delete;

  ValueType getNuclearCoulombEnergy() const { return OneBodyHamiltonian::NuclearCoulombEnergy; }

  boost::multi::array<ValueType, 2> getH1() const { return OneBodyHamiltonian::getH1(); }

  HamiltonianOperations getHamiltonianOperations(bool pureSD,
                                                 bool addCoulomb,
                                                 WALKER_TYPES type,
                                                 std::vector<PsiT_Matrix>& PsiT,
                                                 double cutvn,
                                                 double cutv2,
                                                 TaskGroup_& TGprop,
                                                 TaskGroup_& TGwfn,
                                                 hdf_archive& dump);

  ValueType H(IndexType I, IndexType J) const { return OneBodyHamiltonian::H(I, J); }

  // this should never be used outside initialization routines.
  ValueType H(IndexType I, IndexType J, IndexType K, IndexType L) const
  {
    APP_ABORT("Error: Calling H(I,J,K,L) in THCHamiltonian. \n");
    return ValueType(0.0);
  }

protected:
  // for hamiltonian distribution
  TaskGroup_& TG;

  std::string fileName;

  std::string batched;

  std::string ooc;

  int max_memory_MB;
};

} // namespace afqmc
} // namespace qmcplusplus

#endif
