
#include "QMCDrivers/WFOpt/QMCWFOptFactoryNew.h"
#include "QMCDrivers/WFOpt/QMCOptimize.h"
#include "QMCDrivers/QMCDriverInput.h"
#include "QMCDrivers/VMC/VMCDriverInput.h"
#include "QMCDrivers/WFOpt/QMCOptimizeBatched.h"
#include "QMCDrivers/WFOpt/QMCFixedSampleLinearOptimizeBatched.h"

namespace qmcplusplus
{
QMCOptimizeBatched* QMCWFOptFactoryNew(xmlNodePtr cur,
                                       const int qmc_counter,
                                       MCWalkerConfiguration& w,
                                       TrialWaveFunction& psi,
                                       QMCHamiltonian& h,
                                       HamiltonianPool& hpool,
                                       WaveFunctionPool& wf_pool,
                                       MCPopulation& pop,
                                       SampleStack& samples,
                                       Communicate* comm)
{
  QMCDriverInput qmcdriver_input(qmc_counter);
  qmcdriver_input.readXML(cur);
  VMCDriverInput vmcdriver_input(qmc_counter);
  vmcdriver_input.readXML(cur);

  QMCOptimizeBatched* opt = new QMCOptimizeBatched(w, psi, h, hpool, wf_pool, std::move(qmcdriver_input),
                                                   std::move(vmcdriver_input), pop, samples, comm);
  return opt;
}

QMCFixedSampleLinearOptimizeBatched* QMCWFOptLinearFactoryNew(xmlNodePtr cur,
                                                              const int qmc_counter,
                                                              MCWalkerConfiguration& w,
                                                              TrialWaveFunction& psi,
                                                              QMCHamiltonian& h,
                                                              HamiltonianPool& hpool,
                                                              WaveFunctionPool& wf_pool,
                                                              MCPopulation& pop,
                                                              SampleStack& samples,
                                                              Communicate* comm)
{
  QMCDriverInput qmcdriver_input(qmc_counter);
  qmcdriver_input.readXML(cur);
  VMCDriverInput vmcdriver_input(qmc_counter);
  vmcdriver_input.readXML(cur);

  QMCFixedSampleLinearOptimizeBatched* opt =
      new QMCFixedSampleLinearOptimizeBatched(w, psi, h, hpool, wf_pool, std::move(qmcdriver_input),
                                              std::move(vmcdriver_input), pop, samples, comm);
  return opt;
}

} // namespace qmcplusplus
