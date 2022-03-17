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

namespace qmcplusplus
{
//----------------------------------------------------------------------------
// ProjectData
//----------------------------------------------------------------------------
// constructors and destructors
ProjectData::ProjectData(const std::string& atitle, ProjectData::DriverVersion driver_version)
    : m_title(atitle), m_host("none"), m_date("none"), m_series(0), m_cur(NULL), max_cpu_secs_(360000), driver_version_(driver_version)
{
  myComm  = OHMMS::Controller;
  if (m_title.empty())
    m_title = getDateAndTime("%Y%m%dT%H%M");
}

void ProjectData::setCommunicator(Communicate* c) { myComm = c; }

bool ProjectData::get(std::ostream& os) const
{
  os << "  Project = " << m_title << "\n";
  os << "  date    = " << getDateAndTime("%Y-%m-%d %H:%M:%S %Z\n");
  os << "  host    = " << m_host << "\n";
  return true;
}

bool ProjectData::put(std::istream& is)
{
  std::string t1;
  while (!is.eof())
  {
    if (isdigit(t1[0]))
      m_series = atoi(t1.c_str());
    is >> t1;
    if (t1 == "series")
      is >> m_series;
    else if (t1 == "host")
      is >> m_host;
    else if (t1 == "date")
      is >> m_date;
    else
      m_title = t1;
  }
  reset();
  return true;
}

void ProjectData::advance()
{
  m_series++;
  reset();
}

void ProjectData::rewind()
{
  if (m_series > 0)
    m_series--;
  reset();
}

/**\fn void ProjectData::reset()
 *\brief Construct the root name with m_title and m_series.
 */
void ProjectData::reset()
{
  //int nproc_g = OHMMS::Controller->size();
  int nproc   = myComm->size();
  int nodeid  = myComm->rank();
  int groupid = myComm->getGroupID();
  char fileroot[256], nextroot[256];

  bool no_gtag = (qmc_common.mpi_groups == 1);
  if (no_gtag) //qnproc_g == nproc)
    sprintf(fileroot, "%s.s%03d", m_title.c_str(), m_series);
  else
    sprintf(fileroot, "%s.g%03d.s%03d", m_title.c_str(), groupid, m_series);

  m_projectmain = fileroot;
  //set the communicator name
  myComm->setName(fileroot);
  if (no_gtag)
  {
    if (nproc > 1)
    {
      sprintf(fileroot, ".s%03d.p%03d", m_series, nodeid);
      sprintf(nextroot, ".s%03d.p%03d", m_series + 1, nodeid);
    }
    else
    {
      sprintf(fileroot, ".s%03d", m_series);
      sprintf(nextroot, ".s%03d", m_series + 1);
    }
  }
  else
  {
    if (nproc > 1)
    {
      sprintf(fileroot, ".g%03d.s%03d.p%03d", groupid, m_series, nodeid);
      sprintf(nextroot, ".g%03d.s%03d.p%03d", groupid, m_series + 1, nodeid);
    }
    else
    {
      sprintf(fileroot, ".g%03d.s%03d", groupid, m_series);
      sprintf(nextroot, ".g%03d.s%03d", groupid, m_series + 1);
    }
  }
  m_projectroot = m_title;
  m_projectroot.append(fileroot);
  m_nextroot = m_title;
  m_nextroot.append(nextroot);
  std::stringstream s;
  s << m_series + 1;
  if (m_cur)
    xmlSetProp(m_cur, (const xmlChar*)"series", (const xmlChar*)(s.str().c_str()));
}

bool ProjectData::PreviousRoot(std::string& oldroot) const
{
  oldroot.erase(oldroot.begin(), oldroot.end());
  if (m_series)
  {
    char fileroot[128];
    //int nproc_g = OHMMS::Controller->size();
    int nproc    = myComm->size();
    int nodeid   = myComm->rank();
    int groupid  = myComm->getGroupID();
    bool no_gtag = (qmc_common.mpi_groups == 1);

    if (no_gtag)
    {
      if (nproc > 1)
        sprintf(fileroot, ".s%03d.p%03d", m_series - 1, nodeid);
      else
        sprintf(fileroot, ".s%03d", m_series - 1);
    }
    else
    {
      if (nproc > 1)
        sprintf(fileroot, ".g%03d.s%03d.p%03d", groupid, m_series - 1, nodeid);
      else
        sprintf(fileroot, ".g%03d.s%03d", groupid, m_series - 1);
    }
    oldroot = m_title;
    oldroot.append(fileroot);
    return true;
  }
  else
  {
    return false;
  }
}

bool ProjectData::put(xmlNodePtr cur)
{
  m_cur   = cur;
  m_title = getXMLAttributeValue(cur, "id");
  const std::string series_str(getXMLAttributeValue(cur, "series"));
  if (!series_str.empty())
    m_series = std::stoi(series_str);

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
      m_host = getHostName();
      const XMLNodeString node_string(m_host);
      node_string.setXMLNodeContent(cur);
    }
    if (cname == "date")
    {
      m_date = getDateAndTime();
      const XMLNodeString node_string(m_date);
      node_string.setXMLNodeContent(cur);
    }
    cur = cur->next;
  }
  ///second, add xml nodes, if missing
  if (m_host == "none")
  {
    m_host = getHostName();
    xmlNewChild(m_cur, m_cur->ns, (const xmlChar*)"host", (const xmlChar*)(m_host.c_str()));
  }
  if (m_date == "none")
  {
    m_date = getDateAndTime();
    xmlNewChild(m_cur, m_cur->ns, (const xmlChar*)"date", (const xmlChar*)(m_date.c_str()));
  }
  reset();
  return true;
}

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
