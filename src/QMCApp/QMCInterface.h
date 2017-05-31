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
    
    



/**@file QMCInterface.h
 * @brief Declaration of QMCInterface class
 */
#ifndef QMCPLUSPLUS_INTERFACEAPPLICATIONS_H
#define QMCPLUSPLUS_INTERFACEAPPLICATIONS_H

#include "QMCApp/QMCDriverFactory.h"
#include "QMCApp/QMCAppBase.h"

namespace qmcplusplus
{

/** @ingroup qmcapp
 * @brief Interface application to perform QMC simulations
 *
 * This is a generalized QMC application which can handle multiple ParticleSet,
 * TrialWaveFunction and QMCHamiltonian objects.
 */
class QMCInterface: public QMCDriverFactory,
  public QMCAppBase
{

public:

  ///constructor
  QMCInterface(Communicate* c);

  ///destructor
  ~QMCInterface();

  bool validateXML();
  bool initialize(int myProc, int numProcs);
  bool SetVMC(double dt, int w, int steps, int nblocks);
  bool SetVMCMultiple(double dt, int w, int steps, int nblocks);
  bool SetRQMCMultiple(double dt, int chains, int steps, int nblocks);
  bool process();
  bool execute();
  void SetPtclPos(int id, double* newR);
  void SetPtclPos( std::string set, int id, double* newR);
//    std::vector<double>* GetData();
//		double GetData( std::string estimator, std::string tag);
private:

  /// pointer for xmlNode to be parsed in qmcDriver->process()
  xmlNodePtr runInfoNode;

  ///flag to indicate that a qmc is the first QMC
  bool FirstQMC;

  ///previous configuration file for next qmc node
  std::string PrevConfigFile;

  ///xml mcwalkerset elements for output
  std::vector<xmlNodePtr> m_walkerset;
  ///xml mcwalkerset read-in elements
  std::vector<xmlNodePtr> m_walkerset_in;

  ///execute &lt;qmc/&gt; element
  bool runQMC(xmlNodePtr cur);

  ///add &lt;mcwalkerset/&gt; elements to continue a run
  bool setMCWalkers(xmlXPathContextPtr cur);

  /** add unique ParticleSet, TrialWaveFunction and QMCHamiltonian elements to Pool objects
   */
  void processPWH(xmlNodePtr cur);
};
}
#endif
