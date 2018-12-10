//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



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
namespace qmcplusplus
{

/** @ingroup qmcapp
 * @brief Base class for QMC applications and utilities
 */
class QMCAppBase
{

public:

  ///constructor
  QMCAppBase();

  ///destructor
  ~QMCAppBase();

  /** parse an input file
   * @param infile file to be parsed.
   * @return true if the input file is a valid xml file
   */
  bool parse(const std::string& infile);

  /** save the xml document
   *
   */
  void saveXml();

  /** validate the input file */
  virtual bool validateXML() = 0;

  /** execute the main function */
  virtual bool execute() = 0;

  std::string &getTitle();

protected:

  ///stack of xml document
  std::stack<Libxml2Document*> XmlDocStack;

  ///project description
  ProjectData myProject;

  ///random number controller
  RandomNumberControl myRandomControl;

  ///open a new document
  bool pushDocument(const std::string& infile);
  ///close the current document
  void popDocument();
};
}
#endif
