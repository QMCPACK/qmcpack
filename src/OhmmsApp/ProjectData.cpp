//////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
using namespace std;

#include "OhmmsApp/ProjectData.h"
#include "Message/Communicate.h"
#include "Platforms/sysutil.h"
#if defined(HAVE_LIBBOOST)
#include "boost/date_time/gregorian/gregorian.hpp" 
#endif

namespace OHMMS {

  //----------------------------------------------------------------------------
  // ProjectData
  //----------------------------------------------------------------------------
  // constructors and destructors
  ProjectData::ProjectData(const char* aname): 
    OhmmsElementBase(aname), 
    m_title("asample"),
    m_user("none"), 
    m_host("none"), 
    m_date("none"),
    m_series(0),
    m_cur(0){ 
  }

  bool ProjectData::get(ostream& os) const
  {
    //os << "<Project ID=\""<<m_title << "\" series=\"" << m_series << "/>" << endl;
    os << "  Project = " << m_title << "\n";
#if defined(HAVE_LIBBOOST)
    typedef boost::gregorian::date date_type;
    typedef boost::gregorian::date_facet date_facet;
    typedef boost::date_time::day_clock<date_type> day_clock_type;

    date_type today(day_clock_type::local_day());
    date_facet* facet(new date_facet("%A %B %d, %Y"));
    os.imbue(std::locale(os.getloc(), facet));
    os << "  date    = " << today << "\n";
#endif
    os << "  host    = " << m_host << "\n";
    os << "  user    = " << m_user << "\n";
    return true;
  }

  bool ProjectData::put(istream& is) {

#if defined(ENABLE_GUI)
    // get the data from window
    wxTextCtrl* temp = (wxTextCtrl*)(wxWindow::FindWindowById(ID_PROJECT_TITLE));
    m_title = temp->GetValue();
  
    temp = (wxTextCtrl*)(wxWindow::FindWindowById(ID_PROJECT_SID));
    wxString t = temp->GetValue();
    m_series = atoi(t.c_str());
  
    temp = (wxTextCtrl*)(wxWindow::FindWindowById(ID_DATE));
    wxDateTime now = wxDateTime::Now();
    m_date = now.FormatISODate();
    temp->SetValue(m_date.c_str());
  
    temp = (wxTextCtrl*)(wxWindow::FindWindowById(ID_HOST));
    m_host = wxGetFullHostName();
    temp->SetValue(m_host.c_str());
  
    temp = (wxTextCtrl*)(wxWindow::FindWindowById(ID_USER_ID));
    m_user = wxGetUserId();
    temp->SetValue(m_user.c_str());
#else
    string t1;
    while(!is.eof()) {
      if(isdigit(t1[0])) m_series = atoi(t1.c_str());
      is >> t1; 
      if(t1 == "series") is >> m_series;
      else if(t1 == "user") is >> m_user;
      else if(t1 == "host")   is >> m_host;
      else if(t1 == "date")   is >> m_date;
      else  m_title = t1;
    }
#endif
    reset();
    return true;
  }

  void ProjectData::advance() {
    m_series++; 
    reset();
  }

  void ProjectData::rewind() {
    if(m_series>0) m_series--; 
    reset();
  }

  /**\fn void ProjectData::reset()
   *\brief Construct the root name with m_title and m_series. 
   */
  void ProjectData::reset() {

    char fileroot[128], nextroot[128];
    int nproc = Controller->ncontexts();
    int nodeid = Controller->mycontext(); 
    if(nproc > 1) {
      sprintf(fileroot,".s%03d.p%03d", m_series,nodeid);
      sprintf(nextroot,".s%03d.p%03d", m_series+1,nodeid);
    } else {
      sprintf(fileroot,".s%03d", m_series);
      sprintf(nextroot,".s%03d", m_series+1);
    }
    m_projectroot = m_title;
    m_projectroot.append(fileroot);
    m_nextroot = m_title;
    m_nextroot.append(nextroot);

    std::stringstream s;
    s << m_series+1;
    if(m_cur)
    xmlSetProp(m_cur, (const xmlChar *) "series", (const xmlChar *)(s.str().c_str()));

  }

  bool ProjectData::PreviousRoot(string& oldroot) const {
    oldroot.erase(oldroot.begin(), oldroot.end());
    if(m_series) {
      char fileroot[128];
      int nproc = Controller->ncontexts();
      int nodeid = Controller->mycontext(); 
      if(nproc > 1) {
	sprintf(fileroot,".s%03d.p%03d", m_series-1,nodeid);
      } else {
	sprintf(fileroot,".s%03d", m_series-1);
      }
      oldroot = m_title;
      oldroot.append(fileroot);
      return true;
    } else {
      return false;
    }
  }  

#if defined(HAVE_LIBXML2)

  bool ProjectData::put(xmlNodePtr cur) {

    m_cur = cur;
    xmlDocPtr doc = cur->doc;
    m_title = (const char*)(xmlGetProp(cur, (const xmlChar *) "id"));
    m_series = atoi((const char*)(xmlGetProp(cur, (const xmlChar *) "series")));

    ///first, overwrite the existing xml nodes
    cur = cur->xmlChildrenNode;
    while (cur != NULL) {
      string cname((const char*)(cur->name));
      if(cname == "user") {
        m_user = getUserName(); 
	xmlNodeSetContent(cur,(const xmlChar*)(m_user.c_str()));
      }
      if (cname == "host") {
        m_host = getHostName(); 
	xmlNodeSetContent(cur,(const xmlChar*)(m_host.c_str()));
      }
      if (cname == "date") {
	m_date = getDateAndTime();
	xmlNodeSetContent(cur,(const xmlChar*)(m_date.c_str()));
      }
      cur = cur->next;
    }

    ///second, add xml nodes, if missing
    if(m_host == "none") {
      m_host =  getHostName();
      xmlNewChild(m_cur,m_cur->ns,(const xmlChar*)"host",
		  (const xmlChar*)(m_host.c_str()));
    }
    if(m_date == "none") {
      m_date = getDateAndTime();
      xmlNewChild(m_cur,m_cur->ns,
                  (const xmlChar*)"date",(const xmlChar*)(m_date.c_str()));
    }
    if(m_user == "none") {
      m_user = getUserName();
      xmlNewChild(m_cur,m_cur->ns,
                  (const xmlChar*)"user",(const xmlChar*)(m_user.c_str()));
    }
    reset();

    return true;
  }

#endif
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
