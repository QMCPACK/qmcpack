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
#include "Host/sysutil.h"
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
ProjectData::ProjectData(int my_group_id,
                         int num_groups,
                         const std::string& atitle,
                         ProjectData::DriverVersion driver_version)
    : my_group_id_(my_group_id),
      num_groups_(num_groups),
      title_(atitle),
      host_("none"),
      date_("none"),
      series_(0),
      cur_(NULL),
      max_cpu_secs_(360000),
      driver_version_(driver_version),
      runtime_options_(RuntimeOptions())
{
  if (title_.empty())
    title_ = getDateAndTime("%Y%m%dT%H%M");
}

ProjectData::ProjectData(const std::string& atitle, ProjectData::DriverVersion driver_version)
    : ProjectData(0, 1, atitle, driver_version)
{}

ProjectData::ProjectData(ProjectData::DriverVersion driver_version) : ProjectData(0, 1, "", driver_version) {}

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
  std::array<char, 256> fileroot;

  bool no_gtag = (num_groups_ == 1);
  int file_len{0};
  if (no_gtag)
    file_len = std::snprintf(fileroot.data(), fileroot.size(), "%s.s%03d", title_.c_str(), series_);
  else
    file_len = std::snprintf(fileroot.data(), fileroot.size(), "%s.g%03d.s%03d", title_.c_str(), my_group_id_, series_);

  if (file_len < 0)
    throw std::runtime_error("Error generating project_main");

  project_main_ = std::string(fileroot.data(), file_len);
  std::stringstream s;
  s << series_ + 1;
  if (cur_)
    xmlSetProp(cur_, (const xmlChar*)"series", (const xmlChar*)(s.str().c_str()));
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
