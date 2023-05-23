//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "Configuration.h"
#include "ProgressReportEngine.h"

namespace qmcplusplus
{
void ReportEngine::echo(xmlNodePtr cur, bool recursive)
{
  if (cur == nullptr)
    return;
  app_debug() << R"(<input node=")" << (const char*)(cur->name) << '"';
  xmlAttrPtr att = cur->properties;
  while (att != nullptr)
  {
    app_debug() << "  " << (const char*)(att->name) << R"(=")" << (const char*)(att->children->content) << '"';
    att = att->next;
  }
  app_debug() << R"(/>\n)";
}

bool ReportEngine::DoOutput = false;

} // namespace qmcplusplus
