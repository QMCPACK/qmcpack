
#ifndef AFQMC_TASK_GROUP_H
#define AFQMC_TASK_GROUP_H

#include <vector>
#include <string>
#include <map>
#include <ctime>
#include <sys/time.h>
#include <cstdlib>
#include <ctype.h>
#include <algorithm>
#include <iostream>
#include <ostream>
#include <mpi.h>

#include "AFQMC/config.h"
#include "AFQMC/Memory/utilities.hpp"
#include "mpi3/communicator.hpp"
#include "mpi3/shared_communicator.hpp"

namespace qmcplusplus
{
namespace afqmc
{
using communicator        = boost::mpi3::communicator;
using shared_communicator = boost::mpi3::shared_communicator;

class GlobalTaskGroup
{
public:
  GlobalTaskGroup(communicator& comm)
      : global_(comm),
        //#ifdef ENABLE_CUDA
        // PAMI_DISABLE_IPC issue
        //            node_(comm.split_shared(boost::mpi3::communicator_type::socket,comm.rank())),
        //#else
        node_(comm.split_shared(comm.rank())),
        //#endif
        core_(comm.split(node_.rank(), comm.rank()))
  {
    //    app_log()<<std::endl
    //             <<"**************************************************************\n"
    //             <<" Setting up Global Task Group " <<std::endl;

    // check for consistency
    int dum = node_.size();

    global_.broadcast_n(&dum, 1, 0);
    if (dum != node_.size())
    {
      app_error() << " Error: Inconsistent number of cores in node: " << dum << " " << node_.size() << " "
                  << global_.rank() << std::endl;
      APP_ABORT(" Error in GlobalTaskGroup::GlobalTaskGroup() \n");
    }

    // consistency checks
    dum = core_.size();
    global_.broadcast_n(&dum, 1, 0);
    if (dum != core_.size())
    {
      app_error() << " Error: Inconsistent number of nodes: " << dum << " " << core_.size() << " " << global_.rank()
                  << std::endl;
      APP_ABORT(" Error in GlobalTaskGroup::GlobalTaskGroup() \n");
    }
    dum = core_.rank();
    node_.broadcast_n(&dum, 1, 0);
    if (dum != core_.rank())
    {
      app_error() << " Error: Inconsistent node number: " << dum << " " << core_.rank() << " " << global_.rank()
                  << std::endl;
      APP_ABORT(" Error in GlobalTaskGroup::GlobalTaskGroup() \n");
    }

    //    app_log()<<"**************************************************************" <<std::endl;
  }

  ~GlobalTaskGroup() = default;

  // Not sure if this is necessary, but I want to keep control of the number of communicators
  GlobalTaskGroup(const GlobalTaskGroup& other) = delete;
  GlobalTaskGroup(GlobalTaskGroup&& other)      = delete;
  GlobalTaskGroup& operator=(const GlobalTaskGroup& other) = delete;
  GlobalTaskGroup& operator=(GlobalTaskGroup&& other) = delete;

  int getGlobalRank() const { return global_.rank(); }

  int getGlobalSize() const { return global_.size(); }

  int getTotalNodes() const { return core_.size(); }

  int getTotalCores() const { return node_.size(); }

  int getNodeID() const { return core_.rank(); }

  int getCoreID() const { return node_.rank(); }

  // over full TG using mpi communicator
  void global_barrier() { global_.barrier(); }

  void node_barrier() { node_.barrier(); }

  communicator& Global() { return global_; }
  shared_communicator& Node() { return node_; }
  communicator& Cores() { return core_; }

private:
  // communicators
  communicator& global_;     // global communicator (e.g. MPI_COMM_WORLD or some sub-partition)
  shared_communicator node_; // communicator with all cores in a SHM node
  communicator core_;        // "horizontal" comm, all cores with the same position on a node
};

class TaskGroup_
{
public:
  TaskGroup_(GlobalTaskGroup& gTG, std::string name, int nn, int nc)
      : tgname(name),
        global_(gTG.Global()),
        node_(gTG.Node()),
        core_(gTG.Cores()),
#if defined(ENABLE_CUDA) || defined(ENABLE_HIP)
        local_tg_(node_.split(node_.rank(), node_.rank())),
        tgrp_(),
        tgrp_cores_(),
        tg_heads_(global_.split(0, global_.rank()))
#else
        local_tg_(node_.split(node_.rank() / ((nc < 1) ? (1) : (std::min(nc, node_.size()))), node_.rank())),
        tgrp_(),
        tgrp_cores_(),
        tg_heads_(global_.split(node_.rank() % ((nc < 1) ? (1) : (std::min(nc, node_.size()))), global_.rank()))
#endif
  {
#if defined(ENABLE_CUDA) || defined(ENABLE_HIP)
    if (nc != 1)
    {
      nc = 1;
      app_log() << " WARNING: Code was compiled with ENABLE_CUDA or ENABLE_HIP, setting ncores=1. \n";
    }
#endif
    setup(nn, nc);
  }

  TaskGroup_(TaskGroup_& other, std::string name, int nn, int nc)
      : tgname(name),
        global_(other.Global()),
        node_(other.Node()),
        core_(other.Cores()),
#if defined(ENABLE_CUDA) || defined(ENABLE_HIP)
        local_tg_(node_.split(node_.rank(), node_.rank())),
        tgrp_(),
        tgrp_cores_(),
        tg_heads_(global_.split(0, global_.rank()))
#else
        local_tg_(node_.split(node_.rank() / ((nc < 1) ? (1) : (std::min(nc, node_.size()))), node_.rank())),
        tgrp_(),
        tgrp_cores_(),
        tg_heads_(global_.split(node_.rank() % ((nc < 1) ? (1) : (std::min(nc, node_.size()))), global_.rank()))
#endif
  {
#if defined(ENABLE_CUDA) || defined(ENABLE_HIP)
    if (nc != 1)
    {
      nc = 1;
      app_log() << " WARNING: Code was compiled with ENABLE_CUDA or ENABLE_HIP, setting ncores=1. \n";
    }
#endif
    setup(nn, nc);
  }

  ~TaskGroup_() = default;

  TaskGroup_(const TaskGroup_& other) = delete;
  TaskGroup_(TaskGroup_&& other)
      : global_(other.global_),
        node_(other.node_),
        core_(other.core_),
        local_tg_(other.node_.split(0, other.node_.rank())), // inefficient, but needed to get around lack of
                                                             // default constructor in shared_communicator
        tgrp_(),
        tgrp_cores_(),
        tg_heads_()
  {
    *this = std::move(other);
  }
  TaskGroup_& operator=(const TaskGroup_& other) = delete;
  TaskGroup_& operator                           =(TaskGroup_&& other)
  {
    if (this != &other)
    {
      tgname                          = other.tgname;
      TG_number                       = other.TG_number;
      number_of_TGs                   = other.number_of_TGs;
      nnodes_per_TG                   = other.nnodes_per_TG;
      next_core_circular_node_pattern = other.next_core_circular_node_pattern;
      prev_core_circular_node_pattern = other.prev_core_circular_node_pattern;
      local_tg_                       = std::move(other.local_tg_);
      tgrp_                           = std::move(other.tgrp_);
      tgrp_cores_                     = std::move(other.tgrp_cores_);
      tg_heads_                       = std::move(other.tg_heads_);
#ifdef ENABLE_CUDA
#ifdef BUILD_AFQMC_WITH_NCCL
      nccl_TGcomm_ = std::move(other.nccl_TGcomm_);
      nccl_Stream_ = std::move(other.nccl_Stream_);
#endif
#endif
    }
    return *this;
  }

  int getGlobalRank() const { return global_.rank(); }

  int getGlobalSize() const { return global_.size(); }

  int getTotalNodes() const { return core_.size(); }

  int getTotalCores() const { return node_.size(); }

  int getNodeID() const { return core_.rank(); }

  int getCoreID() const { return node_.rank(); }

  int getLocalTGRank() const { return local_tg_.rank(); }

  int getNCoresPerTG() const { return local_tg_.size(); }

  int getNGroupsPerTG() const { return nnodes_per_TG; }

// THIS NEEDS TO CHANGE IN GPU CODE!
#if defined(ENABLE_CUDA) || defined(ENABLE_HIP)
  int getLocalGroupNumber() const { return getTGRank(); }
#else
  int getLocalGroupNumber() const { return core_.rank() % nnodes_per_TG; }
#endif

  int getTGNumber() const { return TG_number; }

  int getNumberOfTGs() const { return number_of_TGs; }

  int getTGRank() const { return tgrp_.rank(); }

  int getTGSize() const { return tgrp_.size(); }

  // over full TG using mpi communicator
  void global_barrier() { global_.barrier(); }

  // over local node
  void local_barrier() { local_tg_.barrier(); }

  void node_barrier() { node_.barrier(); }

  void barrier() { tgrp_.barrier(); }

  communicator& Global() { return global_; }
  shared_communicator& Node() { return node_; }
  communicator& Cores() { return core_; }
  shared_communicator& TG_local() { return local_tg_; }
  communicator& TG_heads() { return tg_heads_; }
  communicator& TG() { return tgrp_; }
  communicator& TG_Cores() { return tgrp_cores_; }
#ifdef ENABLE_CUDA
#ifdef BUILD_AFQMC_WITH_NCCL
  cudaStream_t& ncclStream() { return nccl_Stream_; }
  ncclComm_t& ncclTG() { return nccl_TGcomm_; }
#endif
#endif

  int next_core() const { return next_core_circular_node_pattern; }
  int prev_core() const { return prev_core_circular_node_pattern; }

private:
  std::string tgname;

  int TG_number;
  int number_of_TGs;
  int nnodes_per_TG;

  int next_core_circular_node_pattern;
  int prev_core_circular_node_pattern;

  // communicators (from Global)
  communicator& global_;      // global communicator (e.g. MPI_COMM_WORLD or some sub-partition)
  shared_communicator& node_; // communicator with all cores in a SHM node
  communicator& core_;        // "horizontal" comm, all cores with the same position on a node
  // communicators (for this TG)
  shared_communicator local_tg_; // all the cores in a node in the same TG
  communicator tgrp_;
  communicator tgrp_cores_;
  communicator tg_heads_;
#ifdef ENABLE_CUDA
#ifdef BUILD_AFQMC_WITH_NCCL
  cudaStream_t nccl_Stream_;
  ncclComm_t nccl_TGcomm_;
#endif
#endif

  void setup(int nn, int nc)
  {
    //    app_log()<<std::endl
    //             <<"**************************************************************\n"
    //             <<" Setting up Task Group: " <<tgname <<std::endl;

    int ncores_per_TG = (nc < 1) ? (1) : (std::min(nc, node_.size()));

    // now setup local_tg_

    if (node_.size() % ncores_per_TG != 0)
    {
      app_error() << "Found " << node_.size() << " cores per node. " << std::endl;
      app_error() << " Error in TaskGroup setup(): Number of cores per node is not divisible by requested number of "
                     "cores in Task Group."
                  << std::endl;
      APP_ABORT(" Error in TaskGroup::TaskGroup() \n");
    }

    if (local_tg_.size() != ncores_per_TG)
    {
      app_log() << nn << " " << nc << std::endl;
      app_error() << "Problems creating local TG: " << local_tg_.size() << " " << ncores_per_TG << std::endl;
      APP_ABORT(" Error in TaskGroup::TaskGroup() \n");
    }

    int ndevices(number_of_devices());

    if (ndevices == 0)
    {
      nnodes_per_TG = (nn < 1) ? (1) : (std::min(nn, core_.size()));

      // assign groups from different nodes to TGs
      if (core_.size() % nnodes_per_TG != 0)
      {
        app_error() << "Found " << core_.size() << " nodes. " << std::endl;
        app_error() << " Error in TaskGroup setup(): Number of nodes is not divisible by requested number of nodes in "
                       "Task Group."
                    << std::endl;
        APP_ABORT(" Error in TaskGroup_::TaskGroup_() \n");
      }

      int myrow, mycol, nrows, ncols, node_in_TG;
      // setup TG grid
      nrows         = node_.size() / ncores_per_TG;
      ncols         = core_.size() / nnodes_per_TG;
      mycol         = core_.rank() / nnodes_per_TG; // simple square asignment
      node_in_TG    = core_.rank() % nnodes_per_TG;
      myrow         = node_.rank() / ncores_per_TG;
      TG_number     = mycol + ncols * myrow;
      number_of_TGs = nrows * ncols;

      // split communicator
      tgrp_ = global_.split(TG_number, global_.rank());

      if (tgrp_.rank() != node_in_TG * ncores_per_TG + local_tg_.rank())
      {
        app_error() << " Error in TG::setup(): Unexpected TG_rank: " << tgrp_.rank() << " " << node_in_TG << " "
                    << local_tg_.rank() << " " << node_in_TG * ncores_per_TG + local_tg_.rank() << std::endl;
        APP_ABORT("Error in TG::setup(): Unexpected TG_rank \n");
      }

      // define ring communication pattern
      // these are the ranks (in TGcomm) of cores with the same core_rank
      // on nodes with id +-1 (with respect to local node).
      next_core_circular_node_pattern = ((node_in_TG + 1) % nnodes_per_TG) * ncores_per_TG + local_tg_.rank();
      if (node_in_TG == 0)
        prev_core_circular_node_pattern = (nnodes_per_TG - 1) * ncores_per_TG + local_tg_.rank();
      else
        prev_core_circular_node_pattern = ((node_in_TG - 1) % nnodes_per_TG) * ncores_per_TG + local_tg_.rank();
    }
    else
    { // ndevices > 0
      // assign groups from the same node to TGs

      nnodes_per_TG = (nn < 1) ? (1) : (std::min(nn, global_.size()));

      if (ncores_per_TG > 1)
      {
        app_error() << " Error in TaskGroup setup(): ncores > 1 incompatible with ndevices>0." << std::endl;
        APP_ABORT(" Error in TaskGroup_::TaskGroup_() \n");
      }
      // limiting to full nodes if nnodes_per_TG > ndevices
      if (global_.size() % nnodes_per_TG != 0)
      {
        app_error() << "Found " << global_.size() << " MPI tasks. " << std::endl;
        app_error() << " Error in TaskGroup setup(): Number of MPI tasks is not divisible by requested number of "
                       "groups in Task Group."
                    << std::endl;
        APP_ABORT(" Error in TaskGroup_::TaskGroup_() \n");
      }

      // MAM: How do I check that distribution of ranks is made along a node first???
      //      Not sure how to check, but print warning if nnodes>1 and ranks are
      //      distributed over nodes first (e.g. adjacent ranks in Global are in different nodes)
      //      This would kill performance!!!
      TG_number     = global_.rank() / nnodes_per_TG;
      number_of_TGs = global_.size() / nnodes_per_TG;

      // split communicator
      tgrp_ = global_.split(TG_number, global_.rank());

      // define ring communication pattern
      next_core_circular_node_pattern = (tgrp_.rank() + 1) % tgrp_.size();
      if (tgrp_.rank() == 0)
        prev_core_circular_node_pattern = tgrp_.size() - 1;
      else
        prev_core_circular_node_pattern = tgrp_.rank() - 1;
    }

#ifdef ENABLE_CUDA
#ifdef BUILD_AFQMC_WITH_NCCL
    qmc_cuda::cuda_check(cudaStreamCreate(&nccl_Stream_), "cudaStreamCreate(&s)");
    {
      ncclUniqueId id;
      if (tgrp_.rank() == 0)
        ncclGetUniqueId(&id);
      MPI_Bcast((void*)&id, sizeof(id), MPI_BYTE, 0, tgrp_.get());
      NCCLCHECK(ncclCommInitRank(&nccl_TGcomm_, tgrp_.size(), id, tgrp_.rank()));
    }
#endif
#endif

    // split communicator
    tgrp_cores_ = tgrp_.split(getLocalTGRank(), getLocalGroupNumber());
  }
};

class TaskGroupHandler
{
public:
  TaskGroupHandler(afqmc::GlobalTaskGroup& gtg_, int nc) : gTG_(gtg_), ncores(nc) {}

  ~TaskGroupHandler() {}

  TaskGroup_& getTG(int nn)
  {
    if (ncores <= 0)
      APP_ABORT(" Error: Calling TaskGroupHandler::getTG() before setting ncores. \n\n\n");
    auto t = TGMap.find(nn);
    if (t == TGMap.end())
    {
      auto p = TGMap.insert(
          std::make_pair(nn, afqmc::TaskGroup_(gTG_, std::string("TaskGroup_") + std::to_string(nn), nn, ncores)));
      if (!p.second)
        APP_ABORT(" Error: Problems creating new TG in TaskGroupHandler::getTG(int). \n");
      return (p.first)->second;
    }
    return t->second;
  }

  void setNCores(int nc) { ncores = nc; }

  afqmc::GlobalTaskGroup& gTG() { return gTG_; }

private:
  afqmc::GlobalTaskGroup& gTG_;

  int ncores;

  std::map<int, afqmc::TaskGroup_> TGMap;
};


} // namespace afqmc

} // namespace qmcplusplus


#endif
