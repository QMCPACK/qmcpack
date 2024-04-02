//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//                    Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_PROJECTDATA_H__
#define QMCPLUSPLUS_PROJECTDATA_H__

#include <unordered_map>
#include "OhmmsData/libxmldefs.h"
#include "Message/Communicate.h"
#include "Utilities/RuntimeOptions.h"

namespace qmcplusplus
{
/**class ProjectData
 *\brief Encapsulate data for a project
 *
 * Default: title_ = getDateAndTime("%Y%m%dT%H%M")
 *
 * \todo This shouldn't contain MPI information it is only used to calculate the
 *       path information and the communication parameters should just be input params to that call.
 *       This reduce the internal state.
 */
class ProjectData
{
public:
  /** Enum for global scope switch of design from legacy driver based to batch driver based.
 *  This effects more than just which drivers are used. Currently it just effects the meaning of
 *  qmc section vmc, dmc, linear name attributes.
 */
  enum class DriverVersion
  {
    LEGACY,
    BATCH
  };

private:
  inline static const std::unordered_map<std::string, DriverVersion> lookup_input_enum_value{{"legacy",
                                                                                              DriverVersion::LEGACY},
                                                                                             {"batch",
                                                                                              DriverVersion::BATCH},
                                                                                             {"batched",
                                                                                              DriverVersion::BATCH}

  };

public:
  /// constructor
  ProjectData(const std::string& atitle = "", DriverVersion de = DriverVersion::LEGACY);

  bool get(std::ostream& os) const;
  bool put(std::istream& is);
  bool put(xmlNodePtr cur);
  void reset();

  ///increment a series number and reset project_root_
  void advance();

  ///roll-back a series number and reset project_root_ by one
  void rewind();

  void setCommunicator(Communicate* c);

  /** 
   * @brief returns the title of the project
   * <project id="det_qmc_short_sdbatch_vmcbatch_mwalkers" series="0">
   * translate to title_ = "det_qmc_short_sdbatch_vmcbatch_mwalkers"
   */
  const std::string& getTitle() const noexcept;

  /** 
   * @brief returns the projectmain of the project, the series id is incremented at every QMC section
   * <project id="det_qmc_short_sdbatch_vmcbatch_mwalkers" series="0">
   * translate to project_main_ = "det_qmc_short_sdbatch_vmcbatch_mwalkers.s000"
   */
  const std::string& currentMainRoot() const noexcept;

  /** 
   * @brief returns the nextroot of the project, the series id is incremented at every QMC section
   * <project id="det_qmc_short_sdbatch_vmcbatch_mwalkers" series="0">
   * translate to project_main_ = "det_qmc_short_sdbatch_vmcbatch_mwalkers.s001"
   */
  const std::string& nextRoot() const noexcept;

  /** 
   * @brief return the root of the previous sequence
   * @param oldroot is composed by the title_ and series_
   */
  bool previousRoot(std::string& oldroot) const;

  int getSeriesIndex() const noexcept;

  int getMaxCPUSeconds() const noexcept;

  DriverVersion getDriverVersion() const noexcept;

  bool isComplex() const noexcept;

  const RuntimeOptions& getRuntimeOptions() const noexcept;

private:
  static DriverVersion lookupDriverVersion(const std::string& enum_value);

  ///title of the project
  std::string title_;

  ///name of the host where the job is running
  std::string host_;

  ///date when the job is executed
  std::string date_;

  ///main root for all the output engines
  std::string project_main_;

  ///processor-dependent root for all the output engines
  std::string project_root_;

  ///root for the next run
  std::string next_root_;

  ///series index
  int series_;

  ///communicator
  Communicate* my_comm_;

  ///the xml node for <Project/>
  xmlNodePtr cur_;

  ///max cpu seconds
  int max_cpu_secs_;

  // The driver version of the project
  DriverVersion driver_version_;

  // tracks runtime options for the project: mixed precision, complex
  RuntimeOptions runtime_options_;
};
} // namespace qmcplusplus

#endif
