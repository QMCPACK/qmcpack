//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//
// File created by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory 
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_AFQMC_SLATERDETOPERATIONS_HPP
#define QMCPLUSPLUS_AFQMC_SLATERDETOPERATIONS_HPP

#include "AFQMC/config.h"
#include "boost/variant.hpp"

#include "AFQMC/SlaterDeterminantOperations/SlaterDetOperations_shared.hpp"

namespace qmcplusplus
{

namespace afqmc
{

class SlaterDetOperations:
        public boost::variant<SlaterDetOperations_shared<ComplexType>>
        //public boost::variant<SlaterDetOperations_shared<ComplexType>, SlaterDetOperations_shared<RealType>>
{

    public:

    SlaterDetOperations(): variant() {
        app_warning()<<(" WARNING: Building SlaterDetOperations with default constructor. \n");
    }

    explicit SlaterDetOperations(SlaterDetOperations_shared<ComplexType>&& other) : variant(std::move(other)) {}
//    explicit SlaterDetOperations(SlaterDetOperations_shared<RealType>&& other) : variant(std::move(other)) {}

    explicit SlaterDetOperations(SlaterDetOperations_shared<ComplexType> const& other) = delete;
//    explicit SlaterDetOperations(SlaterDetOperations_shared<RealType> const& other) = delete;

    SlaterDetOperations(SlaterDetOperations const& other) = delete;
    SlaterDetOperations(SlaterDetOperations && other) = default;

    SlaterDetOperations& operator=(SlaterDetOperations const& other) = delete;
    SlaterDetOperations& operator=(SlaterDetOperations && other) = default;


    // member functions visible outside the variant
    template<class... Args>
    void MixedDensityMatrix(Args&&... args) {
        boost::apply_visitor(
            [&](auto&& a){a.MixedDensityMatrix(std::forward<Args>(args)...);},
            *this
        );
    }

    template<class... Args>
    void MixedDensityMatrix_noHerm(Args&&... args) {
        boost::apply_visitor(
            [&](auto&& a){a.MixedDensityMatrix_noHerm(std::forward<Args>(args)...);},
            *this
        );
    }

    template<class... Args>
    void MixedDensityMatrixForWoodbury(Args&&... args) {
        boost::apply_visitor(
            [&](auto&& a){a.MixedDensityMatrixForWoodbury(std::forward<Args>(args)...);},
            *this
        );
    }

    template<class... Args>
    void MixedDensityMatrixFromConfiguration(Args&&... args) {
        boost::apply_visitor(
            [&](auto&& a){a.MixedDensityMatrixFromConfiguration(std::forward<Args>(args)...);},
            *this
        );
    }

    template<class... Args>
    void Overlap(Args&&... args) {
        boost::apply_visitor(
            [&](auto&& a){a.Overlap(std::forward<Args>(args)...);},
            *this
        );
    }

    template<class... Args>
    void Overlap_noHerm(Args&&... args) {
        boost::apply_visitor(
            [&](auto&& a){a.Overlap_noHerm(std::forward<Args>(args)...);},
            *this
        );
    }

    template<class... Args>
    void OverlapForWoodbury(Args&&... args) {
        boost::apply_visitor(
            [&](auto&& a){a.OverlapForWoodbury(std::forward<Args>(args)...);},
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

    template<class... Args>
    void Orthogonalize(Args&&... args) {
        boost::apply_visitor(
            [&](auto&& a){a.Orthogonalize(std::forward<Args>(args)...);},
            *this
        );
    }

};

}

}

#endif

