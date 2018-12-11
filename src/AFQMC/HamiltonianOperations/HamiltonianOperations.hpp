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

#ifndef QMCPLUSPLUS_AFQMC_HAMILTONIANOPERATIONS_HPP
#define QMCPLUSPLUS_AFQMC_HAMILTONIANOPERATIONS_HPP

#include<fstream>

#include "AFQMC/config.h"
#include "boost/variant.hpp"

#include "AFQMC/HamiltonianOperations/SparseTensor.hpp"
#include "AFQMC/HamiltonianOperations/THCOps.hpp"

namespace qmcplusplus
{

namespace afqmc
{

namespace dummy
{
/*
 * Empty class to avoid need for default constructed HamiltonianOperations.
 * Throws is any visitor is called. 
 */
class dummy_HOps
{
  public:
  dummy_HOps() {};

  boost::multi_array<ComplexType,2> getOneBodyPropagatorMatrix(TaskGroup_& TG, 
                                                               boost::multi_array<ComplexType,1> const& vMF)
  {
    throw std::runtime_error("calling visitor on dummy_HOps object");
    return boost::multi_array<ComplexType,2>{};
  }
  
  template<class Mat, class MatB>
  void energy(Mat&& E, const MatB& Gc, int k, bool addH1=true, bool addEJ=true, bool addEXX=true) 
  {
    throw std::runtime_error("calling visitor on dummy_HOps object");
  }
 
  template<class MatA, class MatB>
  void vHS(const MatA& X, MatB&& v, double a=1., double c=0.) 
  {
    throw std::runtime_error("calling visitor on dummy_HOps object");
  }

  template<class MatA, class MatB>
  void vbias(const MatA& G, MatB&& v, double a=1., double c=0.) 
  {
    throw std::runtime_error("calling visitor on dummy_HOps object");
  }

  template<class... Args>
  void write2hdf(Args&&... args)
  {
    throw std::runtime_error("calling visitor on dummy_HOps object");
  }

  int local_number_of_cholesky_vectors() const
  {
    throw std::runtime_error("calling visitor on dummy_HOps object");
    return 0;
  }

  int global_number_of_cholesky_vectors() const
  {
    throw std::runtime_error("calling visitor on dummy_HOps object");
    return 0;
  }
       
  bool transposed_G_for_vbias() const 
  {
    throw std::runtime_error("calling visitor on dummy_HOps object");
    return false;
  }

  bool transposed_G_for_E() const
  {
    throw std::runtime_error("calling visitor on dummy_HOps object");
    return false;
  }

  bool transposed_vHS() const
  {
    throw std::runtime_error("calling visitor on dummy_HOps object");
    return false;
  }

  bool distribution_over_cholesky_vectors() const
  {
    throw std::runtime_error("calling visitor on dummy_HOps object");
    return false;
  }

};

}


#ifdef QMC_COMPLEX
class HamiltonianOperations: 
        public boost::variant<dummy::dummy_HOps,THCOps<ValueType>,SparseTensor<ComplexType,ComplexType>>
#else
class HamiltonianOperations: 
        public boost::variant<dummy::dummy_HOps,THCOps<ValueType>,
                                  SparseTensor<RealType,RealType>,
                                  SparseTensor<RealType,ComplexType>,
                                  SparseTensor<ComplexType,RealType>,
                                  SparseTensor<ComplexType,ComplexType>
                             >
#endif
{
#ifndef QMC_COMPLEX
    using STRR = SparseTensor<RealType,RealType>;
    using STRC = SparseTensor<RealType,ComplexType>,
    using STCR = SparseTensor<ComplexType,RealType>,
#endif
    using STCC = SparseTensor<ComplexType,ComplexType>;

    public: 

    HamiltonianOperations(): variant() {}
#ifndef QMC_COMPLEX
    explicit HamiltonianOperations(STRR&& other) : variant(std::move(other)) {}
    explicit HamiltonianOperations(STRC&& other) : variant(std::move(other)) {}
    explicit HamiltonianOperations(STCR&& other) : variant(std::move(other)) {}
#endif
    explicit HamiltonianOperations(STCC&& other) : variant(std::move(other)) {}
    explicit HamiltonianOperations(THCOps<ValueType>&& other) : variant(std::move(other)) {}

#ifndef QMC_COMPLEX
    explicit HamiltonianOperations(STRR const& other) = delete;
    explicit HamiltonianOperations(STRC const& other) = delete;
    explicit HamiltonianOperations(STCR const& other) = delete;
#endif
    explicit HamiltonianOperations(STCC const& other) = delete;
    explicit HamiltonianOperations(THCOps<ValueType> const& other) = delete;

    HamiltonianOperations(HamiltonianOperations const& other) = delete;
    HamiltonianOperations(HamiltonianOperations&& other) = default;

    HamiltonianOperations& operator=(HamiltonianOperations const& other) = delete;
    HamiltonianOperations& operator=(HamiltonianOperations&& other) = default;

    template<class... Args>
    boost::multi_array<ComplexType,2> getOneBodyPropagatorMatrix(Args&&... args) {
        return boost::apply_visitor(
            [&](auto&& a){return a.getOneBodyPropagatorMatrix(std::forward<Args>(args)...);},
            *this
        );
    }

    template<class... Args>
    void write2hdf(Args&&... args) {
        boost::apply_visitor(
            [&](auto&& a){a.write2hdf(std::forward<Args>(args)...);},
            *this
        );
    }

    template<class... Args>
    void energy(Args&&... args) {
        boost::apply_visitor(
            [&](auto&& a){a.energy(std::forward<Args>(args)...);},
            *this
        );
    }

    template<class... Args>
    void vHS(Args&&... args) {
        boost::apply_visitor(
            [&](auto&& s){s.vHS(std::forward<Args>(args)...);},
            *this
        );
    }

    template<class... Args>
    void vbias(Args&&... args) {
        boost::apply_visitor(
            [&](auto&& s){s.vbias(std::forward<Args>(args)...);},
            *this
        );
    }

    int local_number_of_cholesky_vectors() const{
        return boost::apply_visitor(
            [&](auto&& a){return a.local_number_of_cholesky_vectors();},
            *this
        );
    }

    int global_number_of_cholesky_vectors() const{
        return boost::apply_visitor(
            [&](auto&& a){return a.global_number_of_cholesky_vectors();},
            *this
        );
    }
    
    bool distribution_over_cholesky_vectors() const {
        return boost::apply_visitor(
            [&](auto&& a){return a.distribution_over_cholesky_vectors();},
            *this
        );
    }


    bool transposed_G_for_vbias() const {
        return boost::apply_visitor(
            [&](auto&& a){return a.transposed_G_for_vbias();},
            *this
        );
    }

    bool transposed_G_for_E() const {
        return boost::apply_visitor(
            [&](auto&& a){return a.transposed_G_for_E();},
            *this
        );
    }

    bool transposed_vHS() const {
        return boost::apply_visitor(
            [&](auto&& a){return a.transposed_vHS();},
            *this
        );
    }
}; 

}

}

#endif
