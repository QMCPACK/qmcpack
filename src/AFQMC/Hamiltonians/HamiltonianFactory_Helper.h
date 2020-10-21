#ifndef QMCPLUSPLUS_AFQMC_HAMILTONIANFACTORY_HELPER_H
#define QMCPLUSPLUS_AFQMC_HAMILTONIANFACTORY_HELPER_H

#include "Configuration.h"
#include "hdf/hdf_archive.h"

#include "AFQMC/config.h"
#include "AFQMC/Hamiltonians/FactorizedSparseHamiltonian.h"

namespace qmcplusplus
{
namespace afqmc
{
FactorizedSparseHamiltonian::shm_csr_matrix read_V2fact(hdf_archive& dump,
                                                        TaskGroup_& TG,
                                                        int nread,
                                                        int NMO,
                                                        int nvecs,
                                                        double cutoff1bar,
                                                        int int_blocks);

}

} // namespace qmcplusplus

#endif
