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
#ifndef OHMMS_QMC_QMCAPPLICATIONBASE_H
#define OHMMS_QMC_QMCAPPLICATIONBASE_H

#include "OhmmsData/OhmmsElementBase.h"
#include "OhmmsApp/ProjectData.h"
#include "OhmmsApp/RandomNumberControl.h"

namespace ohmmsqmc {

  /** Base class for QMC applications and utilities
   */
  class QMCAppBase {

  public:

    ///constructor
    QMCAppBase(int argc, char** argv);

    ///destructor
    ~QMCAppBase();

    /** parse an input file
     * @param infile file to be parsed. 
     * @return true if the input file is a valid xml file
     */
    bool parse(const string& infile);

    /** save the xml document
     *
     */
    void saveXml();

    /** validate the input file */
    virtual bool validateXML() = 0;

    /** execute the main function */
    virtual bool execute() = 0;

  protected:

    xmlDocPtr m_doc;
    xmlNodePtr m_root;

    ///project description
    OHMMS::ProjectData myProject;

    ///random number controller
    OHMMS::RandomNumberControl myRandomControl;
  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
