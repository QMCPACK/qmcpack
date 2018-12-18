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

#include<fstream>

#include "AFQMC/config.h"
#include "boost/variant.hpp"

#include "AFQMC/Propagators/AFQMCSharedPropagator.h"
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
  dummy_Propagator() {};

  template<class WlkSet>
  void Propagate(int steps, WlkSet& wset, RealType& E1, RealType dt, int fix_bias=1) {
    throw std::runtime_error("calling visitor on dummy object");
  } 

  int getNBackProp() { 
    throw std::runtime_error("calling visitor on dummy_Propagator object");
    return 0; 
  }

  bool hybrid_propagation() { 
    throw std::runtime_error("calling visitor on dummy_Propagator object");
    return false;
  }  

};
}

class Propagator: public boost::variant<dummy::dummy_Propagator,AFQMCSharedPropagator,
                                        AFQMCDistributedPropagatorDistCV,
                                        AFQMCDistributedPropagator>
{
    public: 

    Propagator() { 
      APP_ABORT(" Error: Reached default constructor of Propagator. \n");  
    } 
    explicit Propagator(AFQMCSharedPropagator&& other) : variant(std::move(other)) {}
    explicit Propagator(AFQMCSharedPropagator const& other) = delete;

    explicit Propagator(AFQMCDistributedPropagatorDistCV&& other) : variant(std::move(other)) {}
    explicit Propagator(AFQMCDistributedPropagatorDistCV const& other) = delete;

    explicit Propagator(AFQMCDistributedPropagator&& other) : variant(std::move(other)) {}
    explicit Propagator(AFQMCDistributedPropagator const& other) = delete;

    Propagator(Propagator const& other) = delete; 
    Propagator(Propagator&& other) = default; 

    Propagator& operator=(Propagator const& other) = delete; 
    Propagator& operator=(Propagator&& other) = default; 

    int getNBackProp() {
        return boost::apply_visitor(
            [&](auto&& a){return a.getNBackProp();},
            *this
        );
    }

    template<class... Args>
    void Propagate(Args&&... args) {
        boost::apply_visitor(
            [&](auto&& a){a.Propagate(std::forward<Args>(args)...);},
            *this
        );
    }

    bool hybrid_propagation() { 
        return boost::apply_visitor(
            [&](auto&& a){return a.hybrid_propagation();},
            *this
        );
    }

}; 

}

}

#endif
