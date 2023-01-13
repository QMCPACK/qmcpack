
#include "QMCWFOptFactoryNew.h"
#include "QMCDrivers/QMCDriverInput.h"
#include "QMCDrivers/VMC/VMCDriverInput.h"
#include "QMCDrivers/WFOpt/QMCFixedSampleLinearOptimizeBatched.h"
#include "Estimators/EstimatorInputDelegates.h"

namespace qmcplusplus
{
std::unique_ptr<QMCFixedSampleLinearOptimizeBatched> QMCWFOptLinearFactoryNew(
    xmlNodePtr cur,
    const ProjectData& project_data,
    const std::optional<EstimatorManagerInput>& global_emi,
    WalkerConfigurations& wc,
    MCPopulation&& pop,
    SampleStack& samples,
    Communicate* comm)
{
  app_summary() << "\n========================================"
                   "\n  Reading WFOpt driver XML input section"
                   "\n========================================"
                << std::endl;

  QMCDriverInput qmcdriver_input;
  VMCDriverInput vmcdriver_input;
  try
  {
    qmcdriver_input.readXML(cur);
    vmcdriver_input.readXML(cur);
  }
  catch (const std::exception& e)
  {
    throw UniformCommunicateError(e.what());
  }

  auto opt = std::make_unique<QMCFixedSampleLinearOptimizeBatched>(project_data, std::move(qmcdriver_input), global_emi,
                                                                   std::move(vmcdriver_input), wc, std::move(pop),
                                                                   samples, comm);
  return opt;
}

} // namespace qmcplusplus
