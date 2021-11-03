
#include "QMCWFOptFactoryNew.h"
#include "QMCDrivers/QMCDriverInput.h"
#include "QMCDrivers/VMC/VMCDriverInput.h"
#include "QMCDrivers/WFOpt/QMCFixedSampleLinearOptimizeBatched.h"

namespace qmcplusplus
{

QMCFixedSampleLinearOptimizeBatched* QMCWFOptLinearFactoryNew(xmlNodePtr cur,
                                                              const ProjectData& project_data,
                                                              MCWalkerConfiguration& w,
                                                              MCPopulation&& pop,
                                                              SampleStack& samples,
                                                              Communicate* comm)
{
  app_summary() << "\n========================================"
                   "\n  Reading WFOpt driver XML input section"
                   "\n========================================"
                << std::endl;

  QMCDriverInput qmcdriver_input;
  qmcdriver_input.readXML(cur);
  VMCDriverInput vmcdriver_input;
  vmcdriver_input.readXML(cur);

  QMCFixedSampleLinearOptimizeBatched* opt =
      new QMCFixedSampleLinearOptimizeBatched(project_data, w, std::move(qmcdriver_input), std::move(vmcdriver_input),
                                              std::move(pop), samples, comm);
  return opt;
}

} // namespace qmcplusplus
