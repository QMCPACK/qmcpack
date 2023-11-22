//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File refactored from VMC.h
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_QMCDRIVERINTERFACE_H
#define QMCPLUSPLUS_QMCDRIVERINTERFACE_H

#include <string>
#include <vector>
#include <libxml/parser.h>
#include "QMCDrivers/DriverTraits.h"
#include "QMCDrivers/SimpleFixedNodeBranch.h"
#include "QMCDrivers/SFNBranch.h"
#include "Utilities/RandomGenerator.h"

namespace qmcplusplus
{
class QMCHamiltonian;
class TrialWaveFunction;

/** Creates a common base class pointer for QMCDriver and QMCDriverNew
 *  to share.
 *
 *  Once the "unified" driver is done this should be refactored away
 */
class QMCDriverInterface
{
public:
  using BranchEngineType = SimpleFixedNodeBranch;
  using FullPrecRealType              = QMCTraits::FullPrecRealType;
  virtual bool run()                  = 0;
  virtual bool put(xmlNodePtr cur)    = 0;
  virtual void recordBlock(int block) = 0;

  ///return the i-th random generator
  virtual RandomBase<FullPrecRealType>& getRng(int i) = 0;

  virtual void setStatus(const std::string& aname, const std::string& h5name, bool append) = 0;

  virtual void setUpdateMode(bool pbyp)                                 = 0;
  virtual void add_H_and_Psi(QMCHamiltonian* h, TrialWaveFunction* psi) = 0;
  virtual void putWalkers(std::vector<xmlNodePtr>& wset)                = 0;
  virtual void putTraces(xmlNodePtr txml)                               = 0;
  virtual void requestTraces(bool allow_traces)                         = 0;
  virtual void process(xmlNodePtr cur)                                  = 0;
  virtual QMCRunType getRunType()                                       = 0;
  virtual std::string getEngineName()                                   = 0;
  virtual unsigned long getDriverMode()                                 = 0;
  virtual ~QMCDriverInterface() {}

  virtual void setBranchEngine(std::unique_ptr<BranchEngineType>&& be) {}
  virtual std::unique_ptr<BranchEngineType> getBranchEngine() { return nullptr; }

protected:
  virtual const std::string& get_root_name() const = 0;
};

} // namespace qmcplusplus

#endif
