
#ifndef QMCPLUSPLUS_AFQMC_KPFACTORIZEDHAMILTONIAN_H
#define QMCPLUSPLUS_AFQMC_KPFACTORIZEDHAMILTONIAN_H

#include<iostream>
#include<vector>
#include<map>
#include<fstream>

#include "io/hdf_archive.h"

#include "AFQMC/config.h"
#include "AFQMC/Utilities/taskgroup.h"
#include "AFQMC/Numerics/ma_operations.hpp"

#include "AFQMC/Hamiltonians/OneBodyHamiltonian.hpp"

#include "AFQMC/HamiltonianOperations/HamiltonianOperations.hpp"

namespace qmcplusplus
{

namespace afqmc
{

class KPFactorizedHamiltonian: public OneBodyHamiltonian
{

  public:

  using shmSpMatrix = boost::multi::array<SPComplexType,2,shared_allocator<SPComplexType>>;
  using CMatrix = boost::multi::array<ComplexType,2>;

  KPFactorizedHamiltonian(AFQMCInfo const& info, xmlNodePtr cur,
                          boost::multi::array<ComplexType,2>&& h,
                          TaskGroup_& tg_, ValueType nucE=0, ValueType fzcE=0):
                                    OneBodyHamiltonian(info,std::move(h),nucE,fzcE),
                                    TG(tg_),fileName("")
  {

    if( TG.getNumberOfTGs() > 1 )
        APP_ABORT(" Error: Distributed KPFactorizedHamiltonian not yet implemented.\n");

    std::string str("yes");
    ParameterSet m_param;
    m_param.add(cutoff_cholesky,"cutoff_cholesky","double");
    m_param.add(fileName,"filename","std::string");
    m_param.add(nsampleQ,"nsampleQ","int");
    m_param.put(cur);

  }

  ~KPFactorizedHamiltonian() {}

  KPFactorizedHamiltonian(KPFactorizedHamiltonian const& other) = delete;
  KPFactorizedHamiltonian(KPFactorizedHamiltonian && other) = default;
  KPFactorizedHamiltonian& operator=(KPFactorizedHamiltonian const& other) = delete;
  KPFactorizedHamiltonian& operator=(KPFactorizedHamiltonian && other) = default;

  ValueType getNuclearCoulombEnergy() const { return OneBodyHamiltonian::NuclearCoulombEnergy; }

  boost::multi::array<ComplexType,2> getH1() const{ return OneBodyHamiltonian::getH1(); }

  HamiltonianOperations getHamiltonianOperations(bool pureSD, bool addCoulomb, WALKER_TYPES type,
            std::vector<PsiT_Matrix>& PsiT, double cutvn, double cutv2,
            TaskGroup_& TGprop, TaskGroup_& TGwfn, hdf_archive& dump);

  ValueType H(IndexType I, IndexType J) const
  {  return OneBodyHamiltonian::H(I,J); }

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

  double cutoff_cholesky;

  int nsampleQ = -1;

};

}
}

#endif

