//////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License.  See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by:
// Miguel A. Morales, moralessilva2@llnl.gov
//    Lawrence Livermore National Laboratory
//
// File created by:
// Miguel A. Morales, moralessilva2@llnl.gov
//    Lawrence Livermore National Laboratory
////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_AFQMC_PROPAGATOR_HPP
#define QMCPLUSPLUS_AFQMC_PROPAGATOR_HPP

#include <fstream>

#include "AFQMC/config.h"
#include "boost/variant.hpp"

#include "AFQMC/Propagators/AFQMCBasePropagator.h"
#include "AFQMC/Propagators/AFQMCDistributedPropagatorDistCV.h"
#include "AFQMC/Propagators/AFQMCDistributedPropagator.h"

namespace qmcplusplus
{
namespace afqmc
{
namespace dummy
{
/*
 * Empty class to avoid need for default constructed Propagators.
 * Throws is any visitor is called. 
 */
class dummy_Propagator
{
public:
  dummy_Propagator(){};

  template<class WlkSet>
  void Propagate(int steps, WlkSet& wset, RealType& E1, RealType dt, int fix_bias = 1)
  {
    throw std::runtime_error("calling visitor on dummy object");
  }

  template<class WlkSet, class CTens, class CMat>
  void BackPropagate(int steps, int nStabalize, WlkSet& wset, CTens&& Refs, CMat&& detR)
  {
    throw std::runtime_error("calling visitor on dummy object");
  }

  bool hybrid_propagation()
  {
    throw std::runtime_error("calling visitor on dummy_Propagator object");
    return false;
  }

  bool free_propagation()
  {
    throw std::runtime_error("calling visitor on dummy_Propagator object");
    return false;
  }

  int global_number_of_cholesky_vectors() const
  {
    throw std::runtime_error("calling visitor on dummy_Propagator object");
    return 0;
  }

  void generateP1(int, WALKER_TYPES) { throw std::runtime_error("calling visitor on dummy_Propagator object"); }
};
} // namespace dummy

class Propagator : public boost::variant<dummy::dummy_Propagator,
                                         AFQMCBasePropagator,
                                         AFQMCDistributedPropagatorDistCV,
                                         AFQMCDistributedPropagator>
{
public:
  Propagator() { APP_ABORT(" Error: Reached default constructor of Propagator. \n"); }
  explicit Propagator(AFQMCBasePropagator&& other) : variant(std::move(other)) {}
  explicit Propagator(AFQMCBasePropagator const& other) = delete;

  explicit Propagator(AFQMCDistributedPropagatorDistCV&& other) : variant(std::move(other)) {}
  explicit Propagator(AFQMCDistributedPropagatorDistCV const& other) = delete;

  explicit Propagator(AFQMCDistributedPropagator&& other) : variant(std::move(other)) {}
  explicit Propagator(AFQMCDistributedPropagator const& other) = delete;

  Propagator(Propagator const& other) = delete;
  Propagator(Propagator&& other)      = default;

  Propagator& operator=(Propagator const& other) = delete;
  Propagator& operator=(Propagator&& other) = default;

  template<class... Args>
  void Propagate(Args&&... args)
  {
    boost::apply_visitor([&](auto&& a) { a.Propagate(std::forward<Args>(args)...); }, *this);
  }

  template<class... Args>
  void BackPropagate(Args&&... args)
  {
    boost::apply_visitor([&](auto&& a) { a.BackPropagate(std::forward<Args>(args)...); }, *this);
  }

  template<class... Args>
  void generateP1(Args&&... args)
  {
    boost::apply_visitor([&](auto&& a) { a.generateP1(std::forward<Args>(args)...); }, *this);
  }

  bool hybrid_propagation()
  {
    return boost::apply_visitor([&](auto&& a) { return a.hybrid_propagation(); }, *this);
  }

  bool free_propagation()
  {
    return boost::apply_visitor([&](auto&& a) { return a.free_propagation(); }, *this);
  }

  int global_number_of_cholesky_vectors() const
  {
    return boost::apply_visitor([&](auto&& a) { return a.global_number_of_cholesky_vectors(); }, *this);
  }
};

} // namespace afqmc

} // namespace qmcplusplus

#endif
