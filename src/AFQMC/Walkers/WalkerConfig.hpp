#ifndef QMCPLUSPLUS_AFQMC_WALKERCONFIG_HPP
#define QMCPLUSPLUS_AFQMC_WALKERCONFIG_HPP

// basic definitions used by walker set classes
namespace qmcplusplus {
namespace afqmc {

  using wlk_descriptor = std::array<int,4>;
  using wlk_indices = std::array<int,14>;
  enum walker_data { SM, WEIGHT, PHASE, PSEUDO_ELOC_, E1_, EXX_, EJ_, OVLP, PROPAGATORS, HEAD, TAIL, SMN, COS_FAC, WEIGHT_FAC };

}
}

enum LOAD_BALANCE_ALGORITHM { UNDEFINED_LOAD_BALANCE, SIMPLE, ASYNC };
enum BRANCHING_ALGORITHM { UNDEFINED_BRANCHING, PAIR, COMB, MIN_BRANCH, SERIAL_COMB };

#endif
