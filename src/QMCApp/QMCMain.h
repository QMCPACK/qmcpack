//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



/**@file QMCMain.h
 * @brief Declaration of QMCMain class
 */
#ifndef QMCPLUSPLUS_MAINAPPLICATIONS_H
#define QMCPLUSPLUS_MAINAPPLICATIONS_H

#include "QMCApp/QMCDriverFactory.h"
#include "QMCApp/QMCAppBase.h"

namespace qmcplusplus
{

/** @ingroup qmcapp
 * @brief Main application to perform QMC simulations
 *
 * This is a generalized QMC application which can handle multiple ParticleSet,
 * TrialWaveFunction and QMCHamiltonian objects.
 */
class QMCMain: public QMCDriverFactory,
  public QMCAppBase
{

public:

  ///constructor
  QMCMain(Communicate* c);

  ///destructor
  ~QMCMain();

  bool validateXML();
  bool execute();

private:

  ///flag to indicate that a qmc is the first QMC
  bool FirstQMC;

  ///previous configuration file for next qmc node
  std::string PrevConfigFile;

  ///xml mcwalkerset elements for output
  std::vector<xmlNodePtr> m_walkerset;
  ///xml mcwalkerset read-in elements
  std::vector<xmlNodePtr> m_walkerset_in;
  ///traces xml
  xmlNodePtr traces_xml;
  ///qmc sections
  std::vector<std::pair<xmlNodePtr,bool> > m_qmcaction;
  ///pointer to the last node of the main inputfile
  xmlNodePtr lastInputNode;
  ///execute <qmc/> element
  bool runQMC(xmlNodePtr cur);

  ///add <mcwalkerset/> elements to continue a run
  bool setMCWalkers(xmlXPathContextPtr cur);

  /** add unique ParticleSet, TrialWaveFunction and QMCHamiltonian elements to Pool objects
   * @param cur xml node
   * @return false, if contains qmc actions
   */
  bool processPWH(xmlNodePtr cur);

  /** execute loop **/
  void executeLoop(xmlNodePtr cur);
  /** execute qmc
   * @param cur qmc xml node
   * @param noloop if true, this qmc section is not in a loop.
   * @return true, if a section is successfully executed.
   */
  bool executeQMCSection(xmlNodePtr cur, bool noloop=true);
  ///execute <cmc/> element
  bool executeCMCSection(xmlNodePtr cur);
  ///execute <debug/> element
  bool executeDebugSection(xmlNodePtr cur);

  ///execute <postprocess/> element
  void postprocess(xmlNodePtr cur,int qacur);
   
};
}
#endif
