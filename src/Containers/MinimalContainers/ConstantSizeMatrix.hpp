#ifndef CONSTANTSIZEMATRIX_H
#define CONSTANTSIZEMATRIX_H

#include <vector>
#include <sstream>
#include <exception>

#include "OhmmsPETE/OhmmsMatrix.h"

namespace qmcplusplus
{
/** Drop in replacement for OhmmsMatrix as a storage object
 *  backed by PooledMemory, does not support PETE.
 *  rows are contiguous i.e. it's row major like OhmmsMatrix
 *  Constant size removes worry about synchronizing DataSet in Walker.
 *
 *  Reproduces some *creative* operator semantics from OhmmsMatrix,
 *  be careful.
 *
 *  This is intended as a host object, Not for use on accelerators.
 */
template<class T, typename ALLOC = std::allocator<T>>
class ConstantSizeMatrix
{
public:
  ConstantSizeMatrix(size_t m, size_t n, size_t m_max, size_t n_max, T val = T())
      : m_(m), n_(n), m_max_(m_max), n_max_(n_max), capacity_(n_max * m_max), data_(n_max * m_max, val)
  {
    if (n_max <= 0 || n > n_max || m_max <= 0 || m > m_max)
      throw std::runtime_error("Cannot construct ConstantSizeMatrix with an invalid size");
  }

  ConstantSizeMatrix(const ConstantSizeMatrix& csm)
      : m_(csm.m_), n_(csm.n_), m_max_(csm.m_max_), n_max_(csm.n_max_), capacity_(csm.capacity_), data_(csm.capacity_)
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
      throw std::runtime_error("Constant size vector cannot fit matrix to be copied.");
    data_.assign(rhs.begin(), rhs.end());
  }
  void copy(const ConstantSizeMatrix& rhs)
  {
    if (this->data_.capacity() < rhs.size())
      throw std::runtime_error("Constant size vector cannot fit matrix to be copied.");
    if (&rhs == this)
      throw std::runtime_error("ConstantSizeMatrix.copy(ConstantSizeMatrix& rhs) &rhs == this");
    data_.assign(rhs.data_.begin(), rhs.data_.end());
  }

  ConstantSizeMatrix& operator=(const ConstantSizeMatrix& rhs)
  {
    if (this->data_.capacity() < rhs.size())
      throw std::runtime_error("ConstantSizeMatrix cannot take assignment of larger than max size");
    if (&rhs != this)
    {
      if (rhs.n_max_ == n_max_)
        data_.assign(rhs.data_.begin(), rhs.data_.end());
      else
        throw std::runtime_error("ConstnatSizedMatrix assignment for mismatched n_max not yet supported.");
    }
    return *this;
  }

  template<template<typename, class> class RHS,
           typename TP,
           class ALLOCP,
           typename = IsHostSafe<ALLOC>,
           typename = IsHostSafe<ALLOCP>>
  ConstantSizeMatrix<T, ALLOC>& operator=(const RHS<TP, ALLOCP>& rhs)
  {
    if (this->data_.capacity() < rhs.size())
      throw std::runtime_error("ConstantSizeMatrix cannot be assigned a matrix larger than itself.");
    data_.assign(rhs.begin(), rhs.end());
    return *this;
  }

  template<class TP, class ALLOCP, typename = IsHostSafe<ALLOC>, typename = IsHostSafe<ALLOCP>>
  ConstantSizeMatrix<T, ALLOC>& operator=(const Matrix<TP, ALLOCP>& rhs)
  {
    if (this->data_.capacity() < rhs.size())
    {
      std::basic_ostringstream<char> error;
      error << "ConstantSizeMatrix cannot take assignment of larger than max size (" << n_max_ << " by " << m_max_
            << ") \n";
      throw std::runtime_error(error.str());
    }
    if (this->n_max_ < rhs.cols())
    {
      std::basic_ostringstream<char> error;
      error << "ConstantSizeMatrix cannot be assigned a matrix dimension larger than n_max (" << n_max_ << ")\n";
      throw std::runtime_error(error.str());
    }
    n_ = rhs.cols();
    m_ = rhs.rows();
    for (size_t row = 0; row < rhs.rows(); ++row)
      std::copy(rhs[row], rhs[row] + n_, (*this)[row]);
    return *this;
  }

  /** @} */
  // returns a const pointer of i-th row
  const T* operator[](size_t i) const { return data_.data() + i * n_max_; }

  // returns a pointer of i-th row
  T* operator[](size_t i) { return data_.data() + i * n_max_; }

  // returns reference to i-th value in the first row
  template<typename Allocator = ALLOC, typename = IsHostSafe<Allocator>>
  T& operator()(size_t i)
  {
    // As far as I know no code does this and none should
    if (n_ <= i)
      throw std::runtime_error("ConstantSizeMatrix doesn't support access to second row values through operator()");
    return data_[i];
  }

  template<typename Allocator = ALLOC, typename = IsHostSafe<Allocator>>
  T operator()(size_t i) const
  {
    // As far as I know no code does this and none should
    if (n_ <= i)
      throw std::runtime_error("ConstantSizeMatrix doesn't support access to second row values through operator()");
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

  size_t capacity() { return capacity_; }
  size_t n_capacity() { return n_max_; }

  size_t size() const { return n_ * m_; }
  size_t cols() const { return n_; }
  size_t rows() const { return m_; }

  void resize(size_t m, size_t n)
  {
    if (n * m > capacity())
    {
      std::ostringstream error;
      error << "You cannot resize a constant size matrix beyond its initial max dimensions ( " << m << "," << n << " > "
            << m_max_ << "," << n_max_ << " )\n";
      throw std::domain_error(error.str());
    }

    if (n > n_max_)
    {
      n_max_ = n;
      m_max_ = capacity() / n_max_;
    }

    if (m > m_max_)
    {
      m_max_ = m;
      n_max_ = capacity() / m_max_;
    }

    n_ = n;
    m_ = m;
  }

  auto begin() { return data_.begin(); }
  auto end() { return data_.end(); }
  auto begin() const { return data_.begin(); }
  auto end() const { return data_.end(); }

private:
  size_t m_;
  size_t n_;
  size_t m_max_;
  size_t n_max_;
  size_t capacity_;
  std::vector<T, ALLOC> data_;
};

extern template class ConstantSizeMatrix<float>;
extern template class ConstantSizeMatrix<double>;

} // namespace qmcplusplus


#endif /* CONSTANTSIZEMATRIX_H */
