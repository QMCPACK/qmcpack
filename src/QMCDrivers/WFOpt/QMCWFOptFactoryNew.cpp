
#include "QMCWFOptFactoryNew.h"
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
                                       MCPopulation&& pop,
                                       SampleStack& samples,
                                       Communicate* comm)
{
  QMCDriverInput qmcdriver_input(qmc_counter);
  qmcdriver_input.readXML(cur);
  VMCDriverInput vmcdriver_input;
  vmcdriver_input.readXML(cur);

  QMCOptimizeBatched* opt = new QMCOptimizeBatched(w, psi, h, std::move(qmcdriver_input), std::move(vmcdriver_input),
                                                   std::move(pop), samples, comm);
  return opt;
}

QMCFixedSampleLinearOptimizeBatched* QMCWFOptLinearFactoryNew(xmlNodePtr cur,
                                                              const int qmc_counter,
                                                              MCWalkerConfiguration& w,
                                                              TrialWaveFunction& psi,
                                                              QMCHamiltonian& h,
                                                              MCPopulation&& pop,
                                                              SampleStack& samples,
                                                              Communicate* comm)
{
  QMCDriverInput qmcdriver_input(qmc_counter);
  qmcdriver_input.readXML(cur);
  VMCDriverInput vmcdriver_input;
  vmcdriver_input.readXML(cur);

  QMCFixedSampleLinearOptimizeBatched* opt =
      new QMCFixedSampleLinearOptimizeBatched(w, psi, h, std::move(qmcdriver_input), std::move(vmcdriver_input),
                                              std::move(pop), samples, comm);
  return opt;
}

} // namespace qmcplusplus
