
#include "QMCWFOptFactoryNew.h"
#include "QMCDrivers/WFOpt/QMCOptimize.h"
#include "QMCDrivers/QMCDriverInput.h"
#include "QMCDrivers/VMC/VMCDriverInput.h"
#include "QMCDrivers/WFOpt/QMCOptimizeBatched.h"
#include "QMCDrivers/WFOpt/QMCFixedSampleLinearOptimizeBatched.h"

namespace qmcplusplus
{
QMCOptimizeBatched* QMCWFOptFactoryNew(xmlNodePtr cur,
                                       const ProjectData& project_data,
                                       MCWalkerConfiguration& w,
                                       TrialWaveFunction& psi,
                                       QMCHamiltonian& h,
                                       MCPopulation&& pop,
                                       SampleStack& samples,
                                       Communicate* comm)
{
  QMCDriverInput qmcdriver_input;
  qmcdriver_input.readXML(cur);
  VMCDriverInput vmcdriver_input;
  vmcdriver_input.readXML(cur);

  QMCOptimizeBatched* opt = new QMCOptimizeBatched(project_data, w, psi, h, std::move(qmcdriver_input),
                                                   std::move(vmcdriver_input), std::move(pop), samples, comm);
  return opt;
}

QMCFixedSampleLinearOptimizeBatched* QMCWFOptLinearFactoryNew(xmlNodePtr cur,
                                                              const ProjectData& project_data,
                                                              MCWalkerConfiguration& w,
                                                              TrialWaveFunction& psi,
                                                              QMCHamiltonian& h,
                                                              MCPopulation&& pop,
                                                              SampleStack& samples,
                                                              Communicate* comm)
{
  QMCDriverInput qmcdriver_input;
  qmcdriver_input.readXML(cur);
  VMCDriverInput vmcdriver_input;
  vmcdriver_input.readXML(cur);

  QMCFixedSampleLinearOptimizeBatched* opt =
      new QMCFixedSampleLinearOptimizeBatched(project_data, w, psi, h, std::move(qmcdriver_input),
                                              std::move(vmcdriver_input), std::move(pop), samples, comm);
  return opt;
}

} // namespace qmcplusplus
