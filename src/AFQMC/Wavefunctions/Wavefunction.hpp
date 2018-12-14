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

#ifndef QMCPLUSPLUS_AFQMC_WAVEFUNCTION_HPP
#define QMCPLUSPLUS_AFQMC_WAVEFUNCTION_HPP

#include<fstream>

#include "AFQMC/config.h"
#include "boost/variant.hpp"

#include "AFQMC/Wavefunctions/NOMSD.hpp"

namespace qmcplusplus
{

namespace afqmc
{

namespace dummy
{
/*
 * Empty class to avoid need for default constructed Wavefunctions.
 * Throws is any visitor is called. 
 */
class dummy_wavefunction 
{
  private:
  std::vector<ComplexType> ci;
  std::vector<PsiT_Matrix> orbs;
 
  public:
  dummy_wavefunction() {};

  int size_of_G_for_vbias() const {
    throw std::runtime_error("calling visitor on dummy_wavefunction object");  
    return 0; 
  }
  int local_number_of_cholesky_vectors() const{
    throw std::runtime_error("calling visitor on dummy_wavefunction object");  
    return 0; 
  }
  int global_number_of_cholesky_vectors() const{
    throw std::runtime_error("calling visitor on dummy_wavefunction object");  
    return 0; 
  }
  bool distribution_over_cholesky_vectors() const {
    throw std::runtime_error("calling visitor on dummy_wavefunction object");  
    return false; 
  }

/*
  const std::vector<PsiT_Matrix>& getOrbMat() {
    throw std::runtime_error("calling visitor on dummy_wavefunction object");  
    return orbs; 
  } 
  int getOrbSize () { 
    throw std::runtime_error("calling visitor on dummy_wavefunction object");  
    return 0; 
  }
  const std::vector<ComplexType>& getCiCoeff() { 
    throw std::runtime_error("calling visitor on dummy_wavefunction object");  
    return ci; 
  }
*/

  bool transposed_G_for_vbias() const { 
    throw std::runtime_error("calling visitor on dummy_wavefunction object");
    return false; 
  }

  bool transposed_G_for_E() const { 
    throw std::runtime_error("calling visitor on dummy_wavefunction object");
    return false; 
  }

  bool transposed_vHS() const { 
    throw std::runtime_error("calling visitor on dummy_wavefunction object");
    return false; 
  }

  template<class Vec>
  void vMF(Vec&& v) {
    throw std::runtime_error("calling visitor on dummy_wavefunction object");
  }

  template<class MultiArray1D, class task_group>
  boost::multi_array<ComplexType,2> getOneBodyPropagatorMatrix(task_group& tg, MultiArray1D const& vMF) {
    throw std::runtime_error("calling visitor on dummy_wavefunction object");
    return boost::multi_array<ComplexType,2>{};
  }   

  template<class MatG, class MatA>
  void vbias(const MatG& G, MatA&& v, double a=1.0) {  
    throw std::runtime_error("calling visitor on dummy_wavefunction object");  
  }

  template<class MatX, class MatA>
  void vHS(const MatX& X, MatA&& v, double a=1.0) {
    throw std::runtime_error("calling visitor on dummy_wavefunction object");  
  }

  template<class WSet>
  void Energy(WSet& wset) {
    throw std::runtime_error("calling visitor on dummy_wavefunction object");  
  }

  template<class WlkSet, class Mat, class TVec>
  void Energy(const WlkSet& wset, Mat&& E, TVec&& Ov) {
    throw std::runtime_error("calling visitor on dummy_wavefunction object");  
  }

  template<class WlkSet, class MatG>
  void MixedDensityMatrix(const WlkSet& wset, MatG&& G, bool compact=true, bool transpose=false) {
    throw std::runtime_error("calling visitor on dummy_wavefunction object");  
  }

  template<class WlkSet, class MatG, class TVec>
  void MixedDensityMatrix(const WlkSet& wset, MatG&& G, TVec&& Ov, bool compact=true, bool transpose=false) {
    throw std::runtime_error("calling visitor on dummy_wavefunction object");  
  }

  template<class WlkSet, class MatG>
  void MixedDensityMatrix_for_vbias(const WlkSet& wset, MatG&& G) {
    throw std::runtime_error("calling visitor on dummy_wavefunction object");  
  }

  template<class WlkSet, class TVec>
  void Overlap(const WlkSet& wset, TVec&& Ov) {
    throw std::runtime_error("calling visitor on dummy_wavefunction object");  
  }

  template<class WlkSet>
  void Overlap(WlkSet& wset) {
    throw std::runtime_error("calling visitor on dummy_wavefunction object");  
  }

  template<class WlkSet>
  void Orthogonalize(WlkSet& wset, bool impSamp) {
    throw std::runtime_error("calling visitor on dummy_wavefunction object");  
  }

};
}

class Wavefunction: public boost::variant<dummy::dummy_wavefunction,NOMSD>
{
    public: 

    Wavefunction() { 
      APP_ABORT(" Error: Reached default constructor of Wavefunction. \n");  
    } 
    explicit Wavefunction(NOMSD&& other) : variant(std::move(other)) {}
    explicit Wavefunction(NOMSD const& other) = delete;

    Wavefunction(Wavefunction const& other) = delete; 
    Wavefunction(Wavefunction&& other) = default; 

    Wavefunction& operator=(Wavefunction const& other) = delete; 
    Wavefunction& operator=(Wavefunction&& other) = default; 
    
    int size_of_G_for_vbias() const {
        return boost::apply_visitor(
            [&](auto&& a){return a.size_of_G_for_vbias();},
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

    bool distribution_over_cholesky_vectors() const{
        return boost::apply_visitor(
            [&](auto&& a){return a.distribution_over_cholesky_vectors();},
            *this
        );
    }

    bool transposed_G_for_vbias() const{
        return boost::apply_visitor(
            [&](auto&& a){return a.transposed_G_for_vbias();},
            *this
        );
    }

    bool transposed_G_for_E() const{
        return boost::apply_visitor(
            [&](auto&& a){return a.transposed_G_for_E();},
            *this
        );
    }

    bool transposed_vHS() const{
        return boost::apply_visitor(
            [&](auto&& a){return a.transposed_vHS();},
            *this
        );
    }



/*
    std::vector<PsiT_Matrix> const& getOrbMat() const {
        return boost::apply_visitor(
            [&](auto&& a){return a.getOrbMat();},
            *this
        );
    }

    std::vector<ComplexType> const& getCiCoeff() const {
        return boost::apply_visitor(
            [&](auto&& a){return a.getCiCoeff();},
            *this
        );
    }

    int getOrbSize () const {
        return boost::apply_visitor(
            [&](auto&& a){return a.getOrbSize();},
            *this
        );
    }
*/

    template<class... Args>
    void vMF(Args&&... args) {
        boost::apply_visitor(
            [&](auto&& a){a.vMF(std::forward<Args>(args)...);},
            *this
        );
    }

    template<class... Args>
    void vbias(Args&&... args) {
        boost::apply_visitor(
            [&](auto&& a){a.vbias(std::forward<Args>(args)...);},
            *this
        );
    }

    template<class... Args>
    boost::multi_array<ComplexType,2> getOneBodyPropagatorMatrix(Args&&... args) {
        return boost::apply_visitor(
            [&](auto&& a){return a.getOneBodyPropagatorMatrix(std::forward<Args>(args)...);},
            *this
        );
    }

    template<class... Args>
    void vHS(Args&&... args) {
        boost::apply_visitor(
            [&](auto&& a){a.vHS(std::forward<Args>(args)...);},
            *this
        );
    }

    template<class... Args>
    void Energy(Args&&... args) {
        boost::apply_visitor(
            [&](auto&& a){a.Energy(std::forward<Args>(args)...);},
            *this
        );
    }

    template<class... Args>
    void MixedDensityMatrix(Args&&... args) {
        boost::apply_visitor(
            [&](auto&& a){a.MixedDensityMatrix(std::forward<Args>(args)...);},
            *this
        );
    }

    template<class... Args>
    void MixedDensityMatrix_for_vbias(Args&&... args) {
        boost::apply_visitor(
            [&](auto&& a){a.MixedDensityMatrix_for_vbias(std::forward<Args>(args)...);},
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
