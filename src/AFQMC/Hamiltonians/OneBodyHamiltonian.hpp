
#ifndef QMCPLUSPLUS_AFQMC_ONEBODYHAMILTONIAN_H
#define QMCPLUSPLUS_AFQMC_ONEBODYHAMILTONIAN_H

#include<iostream>
#include<vector> 
#include<map> 
#include<fstream>

#include "AFQMC/config.h"
#include "AFQMC/Numerics/SparseMatrixOperations.h"

namespace qmcplusplus
{

namespace afqmc
{

class OneBodyHamiltonian: public AFQMCInfo
{

  typedef std::vector<s2D<ValueType> >::iterator  s2Dit;

  public:
 
  OneBodyHamiltonian(const AFQMCInfo& info, std::vector<s2D<ValueType> >&& h, 
                     ValueType nuc=0, ValueType frz=0): 
                            AFQMCInfo(info), H1(h),
                            NuclearCoulombEnergy(nuc),FrozenCoreEnergy(frz)
  {
    name = "OneBodyHamiltonian";
  }

  OneBodyHamiltonian(OneBodyHamiltonian const& other) = delete;
  OneBodyHamiltonian(OneBodyHamiltonian && other) = default;
  OneBodyHamiltonian& operator=(OneBodyHamiltonian const& other) = delete;
  OneBodyHamiltonian& operator=(OneBodyHamiltonian && other) = default;

  ~OneBodyHamiltonian() {}

  boost::multi_array<ComplexType,2> getH1() const
  {
    boost::multi_array<ComplexType,2> h(extents[NMO][NMO]);
    std::fill_n(h.origin(),h.num_elements(),ComplexType(0)); 
    for(std::vector<s2D<ValueType> >::const_iterator it = H1.begin(); it != H1.end(); it++) {
      h[ std::get<0>(*it) ][ std::get<1>(*it) ] += std::get<2>(*it);
      if(std::get<0>(*it) != std::get<1>(*it))
        h[ std::get<1>(*it) ][ std::get<0>(*it) ] += myconj(std::get<2>(*it));
    }
    return h;
  } 

  // this should never be used outside initialization routines.
  ValueType H(IndexType I, IndexType J) const {
    if( (I>=NMO && J<NMO) || (I<NMO && J>=NMO) ) return ValueType(0);
    I = (I>=NMO)?(I-NMO):(I); 
    J = (J>=NMO)?(J-NMO):(J); 
    if(I <= J) {
      std::vector<s2D<ValueType> >::const_iterator it = std::lower_bound( H1.begin(), H1.end(), std::forward_as_tuple(I,J,static_cast<ValueType>(0.0)),
        [] (const s2D<ValueType>& lhs, const s2D<ValueType>& rhs) 
          {
            return (bool)(std::get<0>(lhs) < std::get<0>(rhs)) ||
              ( !(bool)(std::get<0>(rhs) < std::get<0>(lhs)) &&
              (bool)(std::get<1>(lhs) < std::get<1>(rhs)) );
          }
        );    

      if (it != H1.end() && std::get<0>(*it) == I && std::get<1>(*it) == J) 
        return std::get<2>(*it);  
      else
        return static_cast<ValueType>(0.0);
    } else {
      std::vector<s2D<ValueType> >::const_iterator it = std::lower_bound( H1.begin(), H1.end(), std::forward_as_tuple(J,I,static_cast<ValueType>(0.0)),
        [] (const s2D<ValueType>& lhs, const s2D<ValueType>& rhs) 
          {
            return (bool)(std::get<0>(lhs) < std::get<0>(rhs)) ||
              ( !(bool)(std::get<0>(rhs) < std::get<0>(lhs)) &&
              (bool)(std::get<1>(lhs) < std::get<1>(rhs)) );
          }
        );    

      if (it != H1.end() && std::get<0>(*it) == J && std::get<1>(*it) == I)
        return myconj(std::get<2>(*it));
      else
        return static_cast<ValueType>(0.0);
    }
  }

  ValueType getNuclearCoulombEnergy() const { return NuclearCoulombEnergy; }

  protected:

  // stores one body integrals in s2D format 
  std::vector<s2D<ValueType> > H1;

  // nuclear coulomb term
  ValueType NuclearCoulombEnergy;
  ValueType FrozenCoreEnergy;

};
}
}
#endif

