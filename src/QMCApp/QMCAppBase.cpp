//////////////////////////////////////////////////////////////////
// (c) Copyright 2003- by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim
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
//   Department of Physics, Ohio State University
//   Ohio Supercomputer Center
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#include "QMCApp/QMCAppBase.h"
#include "Utilities/OhmmsInfo.h"

namespace qmcplusplus {

  QMCAppBase::QMCAppBase(int argc, char** argv): 
    m_doc(NULL),m_root(NULL)
  {
  }

  QMCAppBase::~QMCAppBase() {
    if(m_doc != NULL) {
      xmlFreeDoc(m_doc);
    }
    DEBUGMSG("QMCAppBase::~QMCAppBase")
  }

  /** parse an input file
   * @param infile name of an input file
   * @return true, if the document is a valid xml file.
   *
   * The data members m_doc and m_root point to the "current" document and 
   * root element.
   */
  bool QMCAppBase::parse(const string& infile) {

    //clear the context and document
    if(m_doc != NULL) xmlFreeDoc(m_doc);

    // build an XML tree from a the file;
    m_doc = xmlParseFile(infile.c_str());
    if (m_doc == NULL) {
      ERRORMSG("File " << infile << " is invalid")
      return false;
    }

    // Check the document is of the right kind
    xmlNodePtr cur = xmlDocGetRootElement(m_doc);
    if (cur == NULL) {
      ERRORMSG("Empty document");
      return false;
    }

    InFileRoot = string(infile,0,infile.size()-4);

    //set the root and create the context map
    m_root = cur;
    return true;
  }

  void QMCAppBase::saveXml() {
    string newxml(myProject.CurrentRoot());
    //myProject.PreviousRoot(newxml);
    //myProject.rewind();
    newxml.append(".cont.xml");
    LOGMSG("A new xml input file : " << newxml)
    xmlSaveFormatFile(newxml.c_str(),m_doc,1);
  }
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
