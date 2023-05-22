#include "config.h"

namespace qmcplusplus
{
TimerList_t AFQMCTimers;

const std::vector<std::string> AFQMCTimerNames = {"Block",
                                                  "PseudoEnergy",
                                                  "Energy",
                                                  "vHS",
                                                  "X",
                                                  "vbias",
                                                  "G_for_vbias",
                                                  "Propagate",
                                                  "BackPropagate",
                                                  "Energy_comm_overhead",
                                                  "vHS_comm_overhead",
                                                  "population_control",
                                                  "walker_orthogonalization",
                                                  "setup",
                                                  "extra",
                                                  "T1_t",
                                                  "T2_t",
                                                  "T3_t",
                                                  "T4_t",
                                                  "T5_t",
                                                  "T6_t",
                                                  "T7_t",
                                                  "T8_t"};

} // namespace qmcplusplus
