#ifndef QMCPLUSPLUS_AFQMC_WALKERCONFIG_HPP
#define QMCPLUSPLUS_AFQMC_WALKERCONFIG_HPP

// basic definitions used by walker set classes
namespace qmcplusplus
{
namespace afqmc
{
// wlk_descriptor: [ nmo, naea, naeb, nback_prop, nCV, nRefs, nHist]
using wlk_descriptor = std::array<int, 8>;
using wlk_indices    = std::array<int, 17>;
enum walker_data
{
  SM,
  WEIGHT,
  PHASE,
  PSEUDO_ELOC_,
  E1_,
  EXX_,
  EJ_,
  OVLP,
  SMN,
  SM_AUX,
  FIELDS,
  WEIGHT_FAC,
  WEIGHT_HISTORY
};

} // namespace afqmc
} // namespace qmcplusplus

enum LOAD_BALANCE_ALGORITHM
{
  UNDEFINED_LOAD_BALANCE,
  SIMPLE,
  ASYNC
};
enum BRANCHING_ALGORITHM
{
  UNDEFINED_BRANCHING,
  PAIR,
  COMB,
  MIN_BRANCH,
  SERIAL_COMB
};

#endif
