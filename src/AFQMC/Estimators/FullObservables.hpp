i//////////////////////////////////////////////////////////////////////
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

#ifndef QMCPLUSPLUS_AFQMC_FULLOBSERVABLES_HPP
#define QMCPLUSPLUS_AFQMC_FULLOBSERVABLES_HPP

#include "AFQMC/config.h"
#include "boost/variant.hpp"

#include "AFQMC/Estimators/FullObservables_shared.hpp"

namespace qmcplusplus
{

namespace afqmc
{

namespace dummy
{

class dummy_obs
{
  public:
  dummy_obs() {};

  template<class... Args>
  void accumulate(Args&&... args) {
    throw std::runtime_error("calling visitor on dummy_obs object");
  }

  template<class... Args>
  void print(Args&&... args) { 
    throw std::runtime_error("calling visitor on dummy_obs object"); 
  }

};

}

// General interface for observables.
// Given a walker set and an array of (back propagated) slater matrices, 
// this routine will calculate an accumulate all requested observables.
// To make the implementation of the BackPropagated class cleaner, 
// this class also handles all the hdf5 I/O (given a hdf archive).
// This also eliminates the need to move references to arrays between this class
// and the BackPropagated object
class FullObservables: public boost::variant<dummy::dummy_obs,FullObservables_shared> //, FullObservables_batched>
{

  public:

    FullObservables() {
      APP_ABORT(" Error: Reached default constructor of FullObservables().");
    }

    explicit FullObservables(FullObservables_shared&& other) : variant(std::move(other)) {}
    explicit FullObservables(FullObservables_shared const& other) = delete;

    explicit FullObservables(FullObservables_batched&& other) : variant(std::move(other)) {}
    explicit FullObservables(FullObservables_batched const& other) = delete;

    FullObservables(FullObservables const& other) = delete;
    FullObservables(FullObservables&& other) = default;

    FullObservables& operator=(FullObservables const& other) = delete;
    FullObservables& operator=(FullObservables&& other) = default;

    template<class... Args>
    void accumulate(Args&&... args) {
        boost::apply_visitor(
            [&](auto&& a){a.accumulate(std::forward<Args>(args)...);},
            *this
        );
    }

    template<class... Args>
    void print(Args&&... args) {
        boost::apply_visitor(
            [&](auto&& a){a.print(std::forward<Args>(args)...);},
            *this
        );
    }

};


} // afqmc

} // qmcplusplus


#endif
