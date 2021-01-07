#ifndef QMCPLUSPLUS_AFQMC_ESTIMATORBASE_H
#define QMCPLUSPLUS_AFQMC_ESTIMATORBASE_H

#include "AFQMC/config.h"
#include <vector>
#include <iostream>
#include <fstream>

#include "hdf/hdf_multi.h"
#include "hdf/hdf_archive.h"
#include "OhmmsData/libxmldefs.h"

#include "AFQMC/Walkers/WalkerSet.hpp"

namespace qmcplusplus
{
namespace afqmc
{
class EstimatorBase : public AFQMCInfo
{
public:
  EstimatorBase(AFQMCInfo& info) : AFQMCInfo(info) {}

  virtual ~EstimatorBase() {}

  virtual void accumulate_block(WalkerSet& wlks) = 0;

  virtual void accumulate_step(WalkerSet& wlks, std::vector<ComplexType>& curData) = 0;

  virtual void print(std::ofstream& out, hdf_archive& dump, WalkerSet& wlks) = 0;

  virtual void print_timers(std::ofstream& out) {}

  virtual void tags(std::ofstream& out) = 0;

  virtual void tags_timers(std::ofstream& out) {}

  virtual double getEloc() { return 0; }

  virtual double getEloc_step() { return 0; }
};
} // namespace afqmc
} // namespace qmcplusplus

#endif
