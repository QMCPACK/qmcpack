#ifndef QMCPLUSPLUS_OPTIMIZER_TYPE_HEADER
#define QMCPLUSPLUS_OPTIMIZER_TYPE_HEADER
namespace qmcplusplus
{
enum class OptimizerType
{
  LEGACY,
  ONESHIFTONLY,
  ADAPTIVE,
  DESCENT,
  HYBRID
};

const std::map<std::string, OptimizerType> OptimizerNames = {{"Legacy", OptimizerType::LEGACY},
                                                             {"OneShiftOnly", OptimizerType::ONESHIFTONLY},
                                                             {"adaptive", OptimizerType::ADAPTIVE},
                                                             {"descent", OptimizerType::DESCENT},
                                                             {"hybrid", OptimizerType::HYBRID}};

} // namespace qmcplusplus
#endif
