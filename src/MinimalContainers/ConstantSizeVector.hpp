#ifndef CONSTANTSIZEVECTOR_H
#define CONSTANTSIZEVECTOR_H

#include <stdexcept>
#include <vector>

#include "OhmmsPETE/OhmmsVector.h"

namespace qmcplusplus
{

/** Drop in replacement for OhmmsVector as a storage object but not 
 *  for PETE. 
 *
 *  Max size is set and paid at construction it may not change.
 *  This insures compatibility with PooledMemory use in the Walker in a near 0 state design.
 *
 *  Not for use on accelerators.
 */
template<class T,typename ALLOC = std::allocator<T>>
class ConstantSizeVector
{
private:
  size_t n_max_;

public:
  std::vector<T, ALLOC> data_;

  ConstantSizeVector(size_t n, size_t n_max, T val = T()) : n_max_(n_max), data_(n_max, val)
  {
    if (n_max <= 0 || n > n_max )
      throw std::runtime_error("Cannot construct ConstantSizeVector with invalid size");
    data_.resize(n);
  }

  /**@{*/
  /** Methods for assignment or copy of identically sized or smaller
   *  ConstantSizeVector<T, ALLOC>.
   *
   */
  template<class RHS, typename allocator = ALLOC, typename = IsHostSafe<allocator>>
  void copy(const RHS& rhs)
  {
    if (this->data_.capacity() < rhs.size())
      throw std::runtime_error("Constant size vector cannot fit copy.");
    data_.assign(rhs.begin(), rhs.end());
  }
  void copy(const ConstantSizeVector& rhs)
  {
    if (this->data_.capacity() < rhs.size())
      throw std::runtime_error("Constant size vector cannot fit copy.");
    if (&rhs == this)
      throw std::runtime_error("ConstantSizeVector.copy(ConstantSizeVector& rhs) &rhs == this");
    data_.assign(rhs.data_.begin(), rhs.data_.end());
  }
  
  ConstantSizeVector& operator=(const ConstantSizeVector& rhs)
  {
    if (this->data_.capacity() < rhs.size())
      throw std::runtime_error("ConstantSizeVector cannot take assignment of larger than max size");
    if (&rhs == this)
      throw std::runtime_error("ConstantSizeVector operator=(ConstantSizeVector& rhs) &rhs == this");
    data_.assign(rhs.data_.begin(), rhs.data_.end());
    return *this;
  }
  
  template<class TP, class ALLOCP, template<typename TPP = TP, class ALLOCPP=ALLOCP> class RHS, typename std::enable_if<std::is_same<TP, T>::value>, typename = IsHostSafe<ALLOC>>
  ConstantSizeVector<T, ALLOC>& operator=(const RHS<TP, ALLOCP>& rhs) {
//    static_assert(IsHostSafe<ALLOC>);
    if (this->data_.capacity() < rhs.size())
      throw std::runtime_error("ConstantSizeVector cannot take assignment of larger than max size");
    if (&rhs == this)
      throw std::runtime_error("ConstantSizeVector operator=(ConstantSizeVector& rhs) &rhs == this");
    data_.assign(rhs.begin(), rhs.end());
    return *this;
  }

  template<class RHS, typename = std::enable_if_t<!std::is_same<RHS, ConstantSizeVector<T, ALLOC>>::value>>
  ConstantSizeVector<T, ALLOC>& operator=(const RHS& rhs) {
//    static_assert(IsHostSafe<ALLOC>);
    if (this->data_.capacity() < rhs.size())
      throw std::runtime_error("ConstantSizeVector cannot take assignment of larger than max size");
    data_.assign(rhs.begin(), rhs.end());
    return *this;
  }
  /**}@*/

  template<typename allocator = ALLOC, typename = IsHostSafe<allocator>>
  T& operator()(size_t i)
  {
    return data_[i];
  }

  template<typename allocator = ALLOC, typename = IsHostSafe<allocator>>
  const T& operator()(size_t i) const
  {
    return data_[i];
  }

  size_t size() const { return data_.size(); }
  size_t capacity() const { return data_.capacity(); }
  
};

extern template class ConstantSizeVector<double>;

}

#endif /* CONSTANTSIZEVECTOR_H */
