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

#ifndef QMCPLUSPLUS_AFQMC_OBSERVABLE_HPP
#define QMCPLUSPLUS_AFQMC_OBSERVABLE_HPP

#include "AFQMC/config.h"
#include "boost/variant.hpp"

#include "AFQMC/Estimators/Observables/full1rdm.hpp"

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
  void accumulate_block(Args&&... args) {
    throw std::runtime_error("calling visitor on dummy_obs object");
  }

  template<class... Args>
  void accumulate_reference(Args&&... args) {
    throw std::runtime_error("calling visitor on dummy_obs object");
  }

  template<class... Args>
  void print(Args&&... args) { 
    throw std::runtime_error("calling visitor on dummy_obs object"); 
  }

};

}

/*
 * Variant class for observables. 
 * Defines a common interface for all observable classes.
 */
class Observable: public boost::variant<dummy::dummy_obs,full1rdm> 
                                        //,full2rdm,contract1rdm,contract2rdm>
{

  public:

    Observable() {
      APP_ABORT(" Error: Reached default constructor of Observable().");
    }

    explicit Observable(full1rdm && other) : variant(std::move(other)) {}
    explicit Observable(full1rdm const& other) = delete;

/*
    explicit Observable( && other) : variant(std::move(other)) {}
    explicit Observable( const& other) = delete;
*/

    Observable(Observable const& other) = delete;
    Observable(Observable&& other) = default;

    Observable& operator=(Observable const& other) = delete;
    Observable& operator=(Observable&& other) = default;

    template<class... Args>
    void accumulate_block(Args&&... args) {
        boost::apply_visitor(
            [&](auto&& a){a.accumulate_block(std::forward<Args>(args)...);},
            *this
        );
    }

    template<class... Args>
    void accumulate_reference(Args&&... args) {
        boost::apply_visitor(
            [&](auto&& a){a.accumulate_reference(std::forward<Args>(args)...);},
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
