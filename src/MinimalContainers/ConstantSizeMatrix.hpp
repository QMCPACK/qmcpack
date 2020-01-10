#ifndef CONSTANTSIZEMATRIX_H
#define CONSTANTSIZEMATRIX_H

#include "MinimalContainers/ConstantSizeVector.hpp"
namespace qmcplusplus
{
template<class T, typename ALLOC = std::allocator<T>>
class ConstantSizeMatrix
{
public:
  ConstantSizeMatrix(size_t m, size_t n, size_t m_max, size_t n_max, T val = T())
      : m_max_(m_max), n_max_(n_max), data_(n_max * m_max, val)
  {
    if (n_max <= 0 || n > n_max || m_max <= 0 || m > m_max)
      throw std::runtime_error("Cannot construct ConstantSizeMatrix with an invalid size");
  }

  ConstantSizeMatrix(const ConstantSizeMatrix& csm)
      : n_(csm.n_), m_(csm.m_), m_max_(csm.m_max_), n_max_(csm.n_max_), data_(csm.n_max_ * csm.m_max_)
  {
    //I don't just do an = here because of the posible semantics of allocator propagation.
    data_.assign(csm.begin(), csm.end());
  }

    /**@{*/
  /** Methods for assignment or copy of identically sized or smaller
   *  ConstantSizeMatrix<T, ALLOC>.
   *
   */
  template<class RHS, typename allocator = ALLOC, typename = IsHostSafe<allocator>>
  void copy(const RHS& rhs)
  {
    if (this->data_.capacity() < rhs.size())
      throw std::runtime_error("Constant size vector cannot fit copy.");
    data_.assign(rhs.begin(), rhs.end());
  }
  void copy(const ConstantSizeMatrix& rhs)
  {
    if (this->data_.capacity() < rhs.size())
      throw std::runtime_error("Constant size vector cannot fit copy.");
    if (&rhs == this)
      throw std::runtime_error("ConstantSizeMatrix.copy(ConstantSizeMatrix& rhs) &rhs == this");
    data_.assign(rhs.data_.begin(), rhs.data_.end());
  }
  
  ConstantSizeMatrix& operator=(const ConstantSizeMatrix& rhs)
  {
    if (this->data_.capacity() < rhs.size())
      throw std::runtime_error("ConstantSizeMatrix cannot take assignment of larger than max size");
    if (&rhs == this)
      throw std::runtime_error("ConstantSizeMatrix operator=(ConstantSizeMatrix& rhs) &rhs == this");
    data_.assign(rhs.data_.begin(), rhs.data_.end());
    return *this;
  }
  
  template<class TP, class ALLOCP, template<typename TPP = TP, class ALLOCPP=ALLOCP> class RHS, typename std::enable_if<std::is_same<TP, T>::value>, typename = IsHostSafe<ALLOC>>
  ConstantSizeMatrix<T, ALLOC>& operator=(const RHS<TP, ALLOCP>& rhs) {
//    static_assert(IsHostSafe<ALLOC>);
    if (this->data_.capacity() < rhs.size())
      throw std::runtime_error("ConstantSizeMatrix cannot take assignment of larger than max size");
    if (&rhs == this)
      throw std::runtime_error("ConstantSizeMatrix operator=(ConstantSizeMatrix& rhs) &rhs == this");
    data_.assign(rhs.begin(), rhs.end());
    return *this;
  }

  template<class RHS, typename = std::enable_if_t<!std::is_same<RHS, ConstantSizeMatrix<T, ALLOC>>::value>>
  ConstantSizeMatrix<T, ALLOC>& operator=(const RHS& rhs) {
//    static_assert(IsHostSafe<ALLOC>);
    if (this->data_.capacity() < rhs.size())
      throw std::runtime_error("ConstantSizeMatrix cannot take assignment of larger than max size");
    data_.assign(rhs.begin(), rhs.end());
    return *this;
  }
  /**}@*/

  // returns a const pointer of i-th row
  const T* operator[](size_t i) const { return data_.data() + i * n_max_; }

  // returns a pointer of i-th row
  T* operator[](size_t i) { return data_.data() + i * n_max_; }
  
  template<typename Allocator = ALLOC, typename = IsHostSafe<Allocator>>
  T& operator()(size_t i)
  {
    return data_[i];
  }

  // returns the i-th value in n_max_*m_max_ vector
  template<typename Allocator = ALLOC, typename = IsHostSafe<Allocator>>
  T operator()(size_t i) const
  {
    return data_[i];
  }

  // returns val(i,j)
  template<typename Allocator = ALLOC, typename = IsHostSafe<Allocator>>
  T& operator()(size_t i, size_t j)
  {
    return data_[i * n_max_ + j];
  }

  // returns val(i,j)
  template<typename Allocator = ALLOC, typename = IsHostSafe<Allocator>>
  const T& operator()(size_t i, size_t j) const
  {
    return data_[i * n_max_ + j];
  }

  T* data() { return data_.data(); }
  const T* data() const { return data_.data(); }

  size_t sizeofElement() { return sizeof(T); }

  size_t capacity() { return n_max_ * m_max_; }
  size_t n_capacity() { return n_max_; }

  size_t size() const { return n_ * m_; }
  
  void resize(size_t n, size_t m)
  {
    if (n > n_max_ || m > m_max_)
      throw std::runtime_error("You cannot resize a constant size matrix beyond its initial max dimensions");
    n_ = n;
    m_ = m;
  }

  auto begin() { return data_.begin(); }
  auto end() { return data_.end(); }
  auto begin() const { return data_.begin(); }
  auto end() const { return data_.end(); }
private:
  size_t n_;
  size_t m_;
  size_t n_max_;
  size_t m_max_;
  std::vector<T, ALLOC> data_;
};

} // namespace qmcplusplus


#endif /* CONSTANTSIZEMATRIX_H */
