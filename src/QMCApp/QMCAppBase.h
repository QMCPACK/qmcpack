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
#ifndef QMCPLUSPLUS_QMCAPPLICATIONBASE_H
#define QMCPLUSPLUS_QMCAPPLICATIONBASE_H

#include "OhmmsData/OhmmsElementBase.h"
#include "OhmmsData/Libxml2Doc.h"
#include "OhmmsApp/ProjectData.h"
#include "OhmmsApp/RandomNumberControl.h"
#include <stack>
/**@defgroup qmcapp QMC Application Group
 * @brief Application-level classes to manage QMC simulations.
 *
 * The classes in this group are responsble for handling of major xml elements
 * under \<simulation\>.
 */
namespace qmcplusplus {

  /** @ingroup qmcapp
   * @brief Base class for QMC applications and utilities
   */
  class QMCAppBase {

  public:

    ///constructor
    QMCAppBase();

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

    ///stack of xml document
    std::stack<Libxml2Document*> XmlDocStack;

    ///project description
    ProjectData myProject;

    ///random number controller
    RandomNumberControl myRandomControl;

    ///open a new document
    bool pushDocument(const string& infile);
    ///close the current document
    void popDocument();
  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
