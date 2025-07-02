//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License. See LICENSE file in top directory for details.
//
// Copyright (c) 2025 QMCPACK developers.
//
// File developed by: Peter W. Doak, doakpw@ornl.gov, Oak Ridge National
// Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#include "OhmmsData/Libxml2Doc.h"
#include "QMCDrivers/WalkerLogManager.h"
#include "catch.hpp"

#include "Message/Communicate.h"
#include <MockGoldWalkerElements.h>

#include "ValidWalkerLogInput.h"
#include "WalkerLogInput.h"

namespace qmcplusplus
{
class CollectorHolder
{
public:
  void setWalkerLogCollector(UPtr<WalkerLogCollector>&& wlc) { walker_collector_ = std::move(wlc); }
  void startBlock() { walker_collector_->startBlock(); }

private:
  UPtr<WalkerLogCollector> walker_collector_;
};

struct LogAndStuff
{
  WalkerLogManager wlm;
  CollectorHolder ch;
};

TEST_CASE("WalkerLogManager::move", "[drivers]")
{
  Communicate* comm = OHMMS::Controller;
  RuntimeOptions run_time_options;

  auto mgwe     = testing::makeGoldWalkerElementsWithEEEI(comm, run_time_options);
  using WLInput = testing::WalkerLogInputSections;
  Libxml2Document doc;
  bool okay = doc.parseFromString(WLInput::getXml(WLInput::valid::DEFAULT));
  REQUIRE(okay);
  xmlNodePtr node = doc.getRoot();
  WalkerLogInput walker_log_input{node};
  auto make_stuff = [comm](WalkerLogInput& walker_log_input) -> LogAndStuff {
    WalkerLogManager wlm{walker_log_input, true, "root_name", comm};
    CollectorHolder ch;
    ch.setWalkerLogCollector(wlm.makeCollector());
    return {std::move(wlm), std::move(ch)};
  };
  auto l_and_s = make_stuff(walker_log_input);

  // In the past this has resulted in a AddressSanitizer:
  // stack-use-after-return
  // due to access to a dead reference to
  // a state object made back in moved from WalkerLogManager in
  // the factory function.  Seemed to work in release, caught by llvm
  // asan.
  l_and_s.ch.startBlock();
}


} // namespace qmcplusplus
