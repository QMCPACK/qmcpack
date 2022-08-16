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

namespace qmcplusplus
{
// static members
const std::unordered_map<std::string, ProjectData::DriverVersion>
    ProjectData::lookup_input_enum_value_{{"legacy", DriverVersion::LEGACY},
                                          {"batch", DriverVersion::BATCH},
                                          {"batched", DriverVersion::BATCH}};

// PUBLIC

ProjectData::ProjectData(const std::string& atitle, ProjectData::DriverVersion driver_version)
    : title_(atitle),
      host_("none"),
      date_("none"),
      series_index_(0),
      m_cur_(NULL),
      max_cpu_secs_(360000),
      driver_version_(driver_version),
#ifdef QMC_COMPLEX
      is_complex_(true)
#else
      is_complex_(false)
#endif
{
  my_comm_ = OHMMS::Controller;
  if (title_.empty())
    title_ = getDateAndTime("%Y%m%dT%H%M");
}

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
      series_index_ = atoi(t1.c_str());
    is >> t1;
    if (t1 == "series")
      is >> series_index_;
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

bool ProjectData::put(xmlNodePtr cur)
{
  m_cur_ = cur;
  title_ = getXMLAttributeValue(cur, "id");
  const std::string series_str(getXMLAttributeValue(cur, "series"));
  if (!series_str.empty())
    series_index_ = std::stoi(series_str);

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
    xmlNewChild(m_cur_, m_cur_->ns, (const xmlChar*)"host", (const xmlChar*)(host_.c_str()));
  }
  if (date_ == "none")
  {
    date_ = getDateAndTime();
    xmlNewChild(m_cur_, m_cur_->ns, (const xmlChar*)"date", (const xmlChar*)(date_.c_str()));
  }
  reset();
  return true;
}

void ProjectData::reset()
{
  int nproc   = my_comm_->size();
  int nodeid  = my_comm_->rank();
  int groupid = my_comm_->getGroupID();
  char fileroot[256], nextroot[256];

  bool no_gtag = (qmc_common.mpi_groups == 1);

  const int length = no_gtag ? sprintf(fileroot, "%s.s%03d", title_.c_str(), series_index_)
                             : sprintf(fileroot, "%s.g%03d.s%03d", title_.c_str(), groupid, series_index_);

  project_main_ = std::string(fileroot, length);
  //set the communicator name
  my_comm_->setName(fileroot);
  if (no_gtag)
  {
    if (nproc > 1)
    {
      sprintf(fileroot, ".s%03d.p%03d", series_index_, nodeid);
      sprintf(nextroot, ".s%03d.p%03d", series_index_ + 1, nodeid);
    }
    else
    {
      sprintf(fileroot, ".s%03d", series_index_);
      sprintf(nextroot, ".s%03d", series_index_ + 1);
    }
  }
  else
  {
    if (nproc > 1)
    {
      sprintf(fileroot, ".g%03d.s%03d.p%03d", groupid, series_index_, nodeid);
      sprintf(nextroot, ".g%03d.s%03d.p%03d", groupid, series_index_ + 1, nodeid);
    }
    else
    {
      sprintf(fileroot, ".g%03d.s%03d", groupid, series_index_);
      sprintf(nextroot, ".g%03d.s%03d", groupid, series_index_ + 1);
    }
  }
  project_root_ = title_;
  project_root_.append(fileroot);
  next_root_ = title_;
  next_root_.append(nextroot);
  std::stringstream s;
  s << series_index_ + 1;
  if (m_cur_)
    xmlSetProp(m_cur_, (const xmlChar*)"series", (const xmlChar*)(s.str().c_str()));
}

void ProjectData::advance()
{
  series_index_++;
  reset();
}

void ProjectData::rewind()
{
  if (series_index_ > 0)
    series_index_--;
  reset();
}

void ProjectData::setCommunicator(Communicate* c) noexcept { my_comm_ = c; }

const std::string& ProjectData::getTitle() const noexcept { return title_; }

const std::string& ProjectData::currentMainRoot() const noexcept { return project_main_; }

const std::string& ProjectData::nextRoot() const noexcept { return next_root_; }

bool ProjectData::previousRoot(std::string& oldroot) const
{
  oldroot.erase(oldroot.begin(), oldroot.end());
  if (series_index_)
  {
    char fileroot[128];
    const int nproc    = my_comm_->size();
    int nodeid         = my_comm_->rank();
    int groupid        = my_comm_->getGroupID();
    const bool no_gtag = (qmc_common.mpi_groups == 1);

    if (no_gtag)
    {
      if (nproc > 1)
        sprintf(fileroot, ".s%03d.p%03d", series_index_ - 1, nodeid);
      else
        sprintf(fileroot, ".s%03d", series_index_ - 1);
    }
    else
    {
      if (nproc > 1)
        sprintf(fileroot, ".g%03d.s%03d.p%03d", groupid, series_index_ - 1, nodeid);
      else
        sprintf(fileroot, ".g%03d.s%03d", groupid, series_index_ - 1);
    }
    oldroot = title_;
    oldroot.append(fileroot);
    return true;
  }
  else
  {
    return false;
  }
}

int ProjectData::getSeriesIndex() const noexcept { return series_index_; }

int ProjectData::getMaxCPUSeconds() const noexcept { return max_cpu_secs_; }

ProjectData::DriverVersion ProjectData::getDriverVersion() const noexcept { return driver_version_; }

// STATIC PUBLIC
bool ProjectData::isComplex() const noexcept { return is_complex_; }

// PRIVATE
ProjectData::DriverVersion ProjectData::lookupDriverVersion(const std::string& enum_value)
{
  std::string enum_value_str(lowerCase(enum_value));
  try
  {
    return lookup_input_enum_value_.at(enum_value_str);
  }
  catch (std::out_of_range& oor_exc)
  {
    throw UniformCommunicateError("bad_enum_tag_value: " + enum_value_str);
  }
}

} // namespace qmcplusplus
