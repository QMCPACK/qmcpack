
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
    const ParticleSetPool::PoolType& pset_pool,
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

  // This is done so that the application level input structures reflect the actual input to the code.
  // While the actual simulation objects still take singular input structures at construction.
  auto makeEstimatorManagerInput = [](auto& global_emi, auto& local_emi) -> EstimatorManagerInput {
    if (global_emi.has_value() && local_emi.has_value())
      return {global_emi.value(), local_emi.value()};
    else if (global_emi.has_value())
      return {global_emi.value()};
    else if (local_emi.has_value())
      return {local_emi.value()};
    else
      return {};
  };

  // This logic used to be in the QMCFixedSampleLinearOptimizeBatched constructor where it ignored the global_emi because it isn't "really" a driver.
  std::optional<EstimatorManagerInput> dummy_global_emi;
  // In my opinion unwrapping qmcdriver_input_ here and not in QMCDriver new doesn't make sense.
  auto estimator_manager =
      std::make_unique<EstimatorManagerNew>(comm, makeEstimatorManagerInput(dummy_global_emi, qmcdriver_input.get_estimator_manager_input()),
                                            pop.get_golden_hamiltonian(), pop.get_golden_electrons(), pset_pool,
                                            pop.get_golden_twf());


  auto opt =
      std::make_unique<QMCFixedSampleLinearOptimizeBatched>(project_data, std::move(qmcdriver_input),
                                                            std::move(estimator_manager), std::move(vmcdriver_input),
                                                            wc, std::move(pop), samples, comm);
  return opt;
}

} // namespace qmcplusplus
