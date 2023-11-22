//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//                    Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "ProjectData.h"
#include "Message/Communicate.h"
#include "Host/sysutil.h"
#include "Utilities/qmc_common.h"
#include "OhmmsData/ParameterSet.h"
#include "Message/UniformCommunicateError.h"
#include "ModernStringUtils.hpp"

#include <array>

namespace qmcplusplus
{

// PUBLIC
//----------------------------------------------------------------------------
// ProjectData
//----------------------------------------------------------------------------
// constructors and destructors
ProjectData::ProjectData(const std::string& atitle, ProjectData::DriverVersion driver_version)
    : title_(atitle),
      host_("none"),
      date_("none"),
      series_(0),
      cur_(NULL),
      max_cpu_secs_(360000),
      driver_version_(driver_version),
      runtime_options_(RuntimeOptions())
{
  my_comm_ = OHMMS::Controller;
  if (title_.empty())
    title_ = getDateAndTime("%Y%m%dT%H%M");
}

void ProjectData::setCommunicator(Communicate* c) { my_comm_ = c; }

bool ProjectData::get(std::ostream& os) const
{
  os << "  Project = " << title_ << "\n";
  os << "  date    = " << getDateAndTime("%Y-%m-%d %H:%M:%S %Z\n");
  os << "  host    = " << host_ << "\n";
  return true;
}

bool ProjectData::put(std::istream& is)
{
  std::string t1;
  while (!is.eof())
  {
    if (isdigit(t1[0]))
      series_ = atoi(t1.c_str());
    is >> t1;
    if (t1 == "series")
      is >> series_;
    else if (t1 == "host")
      is >> host_;
    else if (t1 == "date")
      is >> date_;
    else
      title_ = t1;
  }
  reset();
  return true;
}

void ProjectData::advance()
{
  series_++;
  reset();
}

void ProjectData::rewind()
{
  if (series_ > 0)
    series_--;
  reset();
}

/**\fn void ProjectData::reset()
 *\brief Construct the root name with title_ and m_series.
 */
void ProjectData::reset()
{
  //int nproc_g = OHMMS::Controller->size();
  int nproc   = my_comm_->size();
  int nodeid  = my_comm_->rank();
  int groupid = my_comm_->getGroupID();
  std::array<char, 256> fileroot;
  std::array<char, 256> nextroot;

  bool no_gtag = (qmc_common.mpi_groups == 1);
  int file_len{0};
  if (no_gtag) //qnproc_g == nproc)
    file_len = std::snprintf(fileroot.data(), fileroot.size(), "%s.s%03d", title_.c_str(), series_);
  else
    file_len = std::snprintf(fileroot.data(), fileroot.size(), "%s.g%03d.s%03d", title_.c_str(), groupid, series_);

  project_main_ = std::string(fileroot.data(), file_len);
  //set the communicator name
  my_comm_->setName(fileroot.data(), file_len);
  int next_len{0};
  if (no_gtag)
  {
    if (nproc > 1)
    {
      file_len = std::snprintf(fileroot.data(), fileroot.size(), ".s%03d.p%03d", series_, nodeid);
      next_len = std::snprintf(nextroot.data(), nextroot.size(), ".s%03d.p%03d", series_ + 1, nodeid);
    }
    else
    {
      file_len = std::snprintf(fileroot.data(), fileroot.size(), ".s%03d", series_);
      next_len = std::snprintf(nextroot.data(), nextroot.size(), ".s%03d", series_ + 1);
    }
  }
  else
  {
    if (nproc > 1)
    {
      file_len = std::snprintf(fileroot.data(), fileroot.size(), ".g%03d.s%03d.p%03d", groupid, series_, nodeid);
      next_len = std::snprintf(nextroot.data(), nextroot.size(), ".g%03d.s%03d.p%03d", groupid, series_ + 1, nodeid);
    }
    else
    {
      file_len = std::snprintf(fileroot.data(), fileroot.size(), ".g%03d.s%03d", groupid, series_);
      next_len = std::snprintf(nextroot.data(), nextroot.size(), ".g%03d.s%03d", groupid, series_ + 1);
    }
  }
  if (file_len < 0)
    throw std::runtime_error("Error generating project_root");
  if (next_len < 0)
    throw std::runtime_error("Error generating next_root");

  project_root_ = title_;
  project_root_.append(fileroot.data(), file_len);
  next_root_ = title_;
  next_root_.append(nextroot.data(), next_len);
  std::stringstream s;
  s << series_ + 1;
  if (cur_)
    xmlSetProp(cur_, (const xmlChar*)"series", (const xmlChar*)(s.str().c_str()));
}

bool ProjectData::previousRoot(std::string& oldroot) const
{
  oldroot.clear();
  if (series_)
  {
    //int nproc_g = OHMMS::Controller->size();
    int nproc    = my_comm_->size();
    int nodeid   = my_comm_->rank();
    int groupid  = my_comm_->getGroupID();
    bool no_gtag = (qmc_common.mpi_groups == 1);
    std::array<char, 128> fileroot;
    int file_len{0};
    if (no_gtag)
    {
      if (nproc > 1)
        file_len = std::snprintf(fileroot.data(), fileroot.size(), ".s%03d.p%03d", series_ - 1, nodeid);
      else
        file_len = std::snprintf(fileroot.data(), fileroot.size(), ".s%03d", series_ - 1);
    }
    else
    {
      if (nproc > 1)
        file_len = std::snprintf(fileroot.data(), fileroot.size(), ".g%03d.s%03d.p%03d", groupid, series_ - 1, nodeid);
      else
        file_len = std::snprintf(fileroot.data(), fileroot.size(), ".g%03d.s%03d", groupid, series_ - 1);
    }
    if (file_len < 0)
      throw std::runtime_error("Error generating olfroot");
    oldroot = title_;
    oldroot.append(fileroot.data(), file_len);
    return true;
  }
  else
  {
    return false;
  }
}

bool ProjectData::put(xmlNodePtr cur)
{
  cur_   = cur;
  title_ = getXMLAttributeValue(cur, "id");
  const std::string series_str(getXMLAttributeValue(cur, "series"));
  if (!series_str.empty())
    series_ = std::stoi(series_str);

  std::string driver_version_str;
  ParameterSet m_param;
  m_param.add(max_cpu_secs_, "max_seconds");
  m_param.add(driver_version_str, "driver_version");
  m_param.put(cur);

  if (!driver_version_str.empty())
    driver_version_ = lookupDriverVersion(driver_version_str);

  ///first, overwrite the existing xml nodes
  cur = cur->xmlChildrenNode;
  while (cur != NULL)
  {
    std::string cname((const char*)(cur->name));
    if (cname == "user")
    {
      // Removed
    }
    if (cname == "host")
    {
      host_ = getHostName();
      const XMLNodeString node_string(host_);
      node_string.setXMLNodeContent(cur);
    }
    if (cname == "date")
    {
      date_ = getDateAndTime();
      const XMLNodeString node_string(date_);
      node_string.setXMLNodeContent(cur);
    }
    cur = cur->next;
  }
  ///second, add xml nodes, if missing
  if (host_ == "none")
  {
    host_ = getHostName();
    xmlNewChild(cur_, cur_->ns, (const xmlChar*)"host", (const xmlChar*)(host_.c_str()));
  }
  if (date_ == "none")
  {
    date_ = getDateAndTime();
    xmlNewChild(cur_, cur_->ns, (const xmlChar*)"date", (const xmlChar*)(date_.c_str()));
  }
  reset();
  return true;
}

const std::string& ProjectData::getTitle() const noexcept { return title_; }

const std::string& ProjectData::currentMainRoot() const noexcept { return project_main_; }

const std::string& ProjectData::nextRoot() const noexcept { return next_root_; }

int ProjectData::getSeriesIndex() const noexcept { return series_; }

int ProjectData::getMaxCPUSeconds() const noexcept { return max_cpu_secs_; }

ProjectData::DriverVersion ProjectData::getDriverVersion() const noexcept { return driver_version_; }

bool ProjectData::isComplex() const noexcept { return runtime_options_.is_complex; }

const RuntimeOptions& ProjectData::getRuntimeOptions() const noexcept { return runtime_options_; }

// PRIVATE
ProjectData::DriverVersion ProjectData::lookupDriverVersion(const std::string& enum_value)
{
  std::string enum_value_str(lowerCase(enum_value));
  try
  {
    return lookup_input_enum_value.at(enum_value_str);
  }
  catch (std::out_of_range& oor_exc)
  {
    throw UniformCommunicateError("bad_enum_tag_value: " + enum_value_str);
  }
}

} // namespace qmcplusplus
