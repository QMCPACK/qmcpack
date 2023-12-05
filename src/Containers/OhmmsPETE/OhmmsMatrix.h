//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////

#ifndef OHMMS_PETE_MATRIX_H
#define OHMMS_PETE_MATRIX_H

#include <cstdlib>
#include <type_traits>
#include <iostream>
#include "PETE/PETE.h"
#include "OhmmsVector.h"

namespace qmcplusplus
{
template<class T, typename Alloc = std::allocator<T>>
class Matrix
{
public:
  using Type_t        = T;
  using value_type    = T;
  using pointer       = T*;
  using const_pointer = const T*;
  using Container_t   = Vector<T, Alloc>;
  using size_type     = typename Container_t::size_type;
  using iterator      = typename Container_t::iterator;
  using This_t        = Matrix<T, Alloc>;

  Matrix() : D1(0), D2(0), TotSize(0) {} // Default Constructor initializes to zero.

  Matrix(size_type n)
  {
    resize(n, n);
    //assign(*this, T());
  }

  Matrix(size_type n, size_type m)
  {
    resize(n, m);
    //assign(*this, T());
  }

  /** constructor with an initialized ref */
  inline Matrix(T* ref, size_type n, size_type m) : D1(n), D2(m), TotSize(n * m), X(ref, n * m) {}

  /** This allows construction of a Matrix on another containers owned memory that is using a dualspace allocator.
   *  It can be any span of that memory.
   *  You're going to get a bunch of compile errors if the Container in questions is not using a the QMCPACK
   *  realspace dualspace allocator "interface"
   */
  template<typename CONTAINER>
  Matrix(CONTAINER& other, T* ref, size_type n, size_type m) : D1(n), D2(m), TotSize(n * m), X(other, ref, n * m)
  {}

  // Copy Constructor
  Matrix(const This_t& rhs)
  {
    resize(rhs.D1, rhs.D2);
    if (qmc_allocator_traits<Alloc>::is_host_accessible)
      assign(*this, rhs);
  }

  // Destructor
  ~Matrix() {}

  inline size_type size() const { return TotSize; }
  inline size_type rows() const { return D1; }
  inline size_type cols() const { return D2; }
  inline size_type size1() const { return D1; }
  inline size_type size2() const { return D2; }
  inline size_type size(int i) const { return (i == 0) ? D1 : D2; }
  inline size_type extent(int i) const { return (i == 0) ? D1 : D2; }

  //   inline const T* begin() const { return X.begin();}
  //   inline T* begin() { return X.begin();}
  //   inline const T* end() const { return X.end();}
  //   inline T* end()   { return X.end();}

  inline typename Container_t::iterator begin() { return X.begin(); }
  inline typename Container_t::iterator end() { return X.end(); }
  inline typename Container_t::const_iterator begin() const { return X.begin(); }
  inline typename Container_t::const_iterator end() const { return X.end(); }

  inline typename Container_t::iterator begin(int i) { return X.begin() + i * D2; }
  inline typename Container_t::const_iterator begin(int i) const { return X.begin() + i * D2; }

  /// Resize the container. For performance consideration, previous data may or may not get kept.
  /// Please avoid relying on previous data after resizing.
  inline void resize(size_type n, size_type m)
  {
    static_assert(std::is_same<value_type, typename Alloc::value_type>::value,
                  "Matrix and Alloc data types must agree!");
    D1      = n;
    D2      = m;
    TotSize = n * m;
    X.resize(n * m);
  }

  // free the matrix storage
  inline void free() { X.free(); }

  // Attach to pre-allocated memory
  inline void attachReference(T* ref) { X.attachReference(ref, TotSize); }

  inline void attachReference(T* ref, size_type n, size_type m)
  {
    D1      = n;
    D2      = m;
    TotSize = n * m;
    X.attachReference(ref, TotSize);
  }

  /** Attach to pre-allocated memory and propagate the allocator of the owning container.
   *  Required for sane access to dual space memory
   */
  template<typename CONTAINER>
  inline void attachReference(const CONTAINER& other, T* ref, size_type n, size_type m)
  {
    D1      = n;
    D2      = m;
    TotSize = n * m;
    X.attachReference(other, ref, TotSize);
  }

  template<typename Allocator = Alloc, typename = IsHostSafe<Allocator>>
  inline void add(size_type n) // you can add rows: adding columns are forbidden
  {
    X.insert(X.end(), n * D2, T());
    D1 += n;
  }

  template<typename Allocator = Alloc, typename = IsHostSafe<Allocator>>
  inline void copy(const This_t& rhs)
  {
    resize(rhs.D1, rhs.D2);
    assign(*this, rhs);
  }

  /** This assigns from a matrix with larger row size (used for alignment)
   *  to whatever the rowsize is here.
   *  Hacky but so is just making the matrix n x (n + padding) to handle row alignment.
   *  This is unoptimized.
   */
  template<class T_FROM, typename ALLOC_FROM>
  void assignUpperLeft(const Matrix<T_FROM, ALLOC_FROM>& from)
  {
    auto& this_ref    = *this;
    const size_t cols = std::min(this_ref.cols(), from.cols());
    const size_t rows = std::min(this_ref.rows(), from.rows());
    for (int i = 0; i < rows; ++i)
      for (int j = 0; j < cols; ++j)
        this_ref(i, j) = from(i, j);
  }

  // Assignment Operators
  inline This_t& operator=(const This_t& rhs)
  {
    resize(rhs.D1, rhs.D2);
    if (qmc_allocator_traits<Alloc>::is_host_accessible)
      assign(*this, rhs);
    return *this;
  }

  template<class RHS, typename Allocator = Alloc, typename = IsHostSafe<Allocator>>
  This_t& operator=(const RHS& rhs)
  {
    return assign(*this, rhs);
  }

  // Get and Set Operations for assignment operators
  // returns a pointer of i-th row
  inline pointer data() { return X.data(); }

  // returns a pointer of i-th row
  inline const_pointer data() const { return X.data(); }

  template<typename Allocator = Alloc, typename = IsDualSpace<Allocator>>
  inline pointer device_data()
  {
    return X.device_data();
  }
  template<typename Allocator = Alloc, typename = IsDualSpace<Allocator>>
  inline const_pointer device_data() const
  {
    return X.device_data();
  }

  // returns a const pointer of i-th row
  inline const Type_t* data(size_type i) const { return X.data() + i * D2; }

  /// returns a pointer of i-th row, g++ iterator problem
  inline Type_t* data(size_type i) { return X.data() + i * D2; }

  inline pointer first_address() { return X.data(); }

  // returns a pointer of i-th row
  inline const_pointer first_address() const { return X.data(); }

  inline pointer last_address() { return X.data() + TotSize; }

  // returns a pointer of i-th row
  inline const Type_t* last_address() const { return X.data() + TotSize; }


  // returns a const pointer of i-th row
  inline const Type_t* operator[](size_type i) const { return X.data() + i * D2; }

  /// returns a pointer of i-th row, g++ iterator problem
  inline Type_t* operator[](size_type i) { return X.data() + i * D2; }

  template<typename Allocator = Alloc, typename = IsHostSafe<Allocator>>
  inline Type_t& operator()(size_type i)
  {
    return X[i];
  }

  // returns the i-th value in D1*D2 vector
  template<typename Allocator = Alloc, typename = IsHostSafe<Allocator>>
  inline Type_t operator()(size_type i) const
  {
    return X[i];
  }

  // returns val(i,j)
  template<typename Allocator = Alloc, typename = IsHostSafe<Allocator>>
  inline Type_t& operator()(size_type i, size_type j)
  {
    return X[i * D2 + j];
  }

  // returns val(i,j)
  template<typename Allocator = Alloc, typename = IsHostSafe<Allocator>>
  inline const Type_t& operator()(size_type i, size_type j) const
  {
    return X[i * D2 + j];
  }

  template<typename Allocator = Alloc, typename = IsHostSafe<Allocator>>
  inline void swap_rows(int r1, int r2)
  {
    for (int col = 0; col < D2; col++)
    {
      Type_t tmp       = (*this)(r1, col);
      (*this)(r1, col) = (*this)(r2, col);
      (*this)(r2, col) = tmp;
    }
  }

  template<typename Allocator = Alloc, typename = IsHostSafe<Allocator>>
  inline void swap_cols(int c1, int c2)
  {
    for (int row = 0; row < D1; row++)
    {
      Type_t tmp       = (*this)(row, c1);
      (*this)(row, c1) = (*this)(row, c2);
      (*this)(row, c2) = tmp;
    }
  }

  template<class IT, typename Allocator = Alloc, typename = IsHostSafe<Allocator>>
  inline void replaceRow(IT first, size_type i)
  {
    std::copy(first, first + D2, X.begin() + i * D2);
  }

  template<class IT, typename Allocator = Alloc, typename = IsHostSafe<Allocator>>
  inline void replaceColumn(IT first, size_type j)
  {
    typename Container_t::iterator ii(X.begin() + j);
    for (int i = 0; i < D1; i++, ii += D2)
      *ii = *first++;
  }

  template<class IT, typename Allocator = Alloc, typename = IsHostSafe<Allocator>>
  inline void add2Column(IT first, size_type j)
  {
    typename Container_t::iterator ii(X.begin() + j);
    for (int i = 0; i < D1; i++, ii += D2)
      *ii += *first++;
  }

  /**
   * \param sub an input array to be copied
   * \param d1  row-dimension of the input array
   * \param d2  column-dimension of the input array
   * \param i0  starting row where the copying is done
   * \param j0  starting column where the copying is done
   */
  template<class T1, typename Allocator = Alloc, typename = IsHostSafe<Allocator>>
  inline void add(const T1* sub, size_type d1, size_type d2, size_type i0, size_type j0)
  {
    int ii = 0;
    for (int i = 0; i < d1; i++)
    {
      int kk = (i0 + i) * D2 + j0;
      for (int j = 0; j < d2; j++)
      {
        X[kk++] += sub[ii++];
      }
    }
  }

  template<class T1, typename Allocator = Alloc, typename = IsHostSafe<Allocator>>
  inline void add(const T1* sub, size_type d1, size_type d2, size_type i0, size_type j0, const T& phi)
  {
    size_type ii = 0;
    for (size_type i = 0; i < d1; i++)
    {
      int kk = (i0 + i) * D2 + j0;
      for (size_type j = 0; j < d2; j++)
      {
        X[kk++] += phi * sub[ii++];
      }
    }
  }

  template<class SubMat_t, typename Allocator = Alloc, typename = IsHostSafe<Allocator>>
  inline void add(const SubMat_t& sub, unsigned int i0, unsigned int j0)
  {
    size_type ii = 0;
    for (size_type i = 0; i < sub.rows(); i++)
    {
      int kk = (i0 + i) * D2 + j0;
      for (size_type j = 0; j < sub.cols(); j++)
      {
        X[kk++] += sub(ii++);
      }
    }
  }

  template<typename Allocator = Alloc, typename = IsHostSafe<Allocator>>
  inline void add(const This_t& sub, unsigned int i0, unsigned int j0)
  {
    size_type ii = 0;
    for (size_type i = 0; i < sub.rows(); i++)
    {
      int kk = (i0 + i) * D2 + j0;
      for (size_type j = 0; j < sub.cols(); j++)
      {
        X[kk++] += sub[ii++];
      }
    }
  }

  template<class Msg, typename Allocator = Alloc, typename = IsHostSafe<Allocator>>
  inline Msg& putMessage(Msg& m)
  {
    m.Pack(X.data(), D1 * D2);
    return m;
  }

  template<class Msg, typename Allocator = Alloc, typename = IsHostSafe<Allocator>>
  inline Msg& getMessage(Msg& m)
  {
    m.Unpack(X.data(), D1 * D2);
    return m;
  }

  // Abstract Dual Space Transfers
  template<typename Allocator = Alloc, typename = IsDualSpace<Allocator>>
  void updateTo()
  {
    X.updateTo();
  }
  template<typename Allocator = Alloc, typename = IsDualSpace<Allocator>>
  void updateFrom()
  {
    X.updateFrom();
  }

protected:
  size_type D1, D2;
  size_type TotSize;
  Container_t X;
};

template<class T, class Alloc>
bool operator==(const Matrix<T, Alloc>& lhs, const Matrix<T, Alloc>& rhs)
{
  static_assert(qmc_allocator_traits<Alloc>::is_host_accessible, "operator== requires host accessible Vector.");
  if (lhs.size() == rhs.size())
  {
    for (int i = 0; i < rhs.size(); i++)
      if (lhs(i) != rhs(i))
        return false;
    return true;
  }
  else
    return false;
}

template<class T, class Alloc>
bool operator!=(const Matrix<T, Alloc>& lhs, const Matrix<T, Alloc>& rhs)
{
  static_assert(qmc_allocator_traits<Alloc>::is_host_accessible, "operator== requires host accessible Vector.");
  return !(lhs == rhs);
}


// I/O
template<class T, typename Alloc>
std::ostream& operator<<(std::ostream& out, const Matrix<T, Alloc>& rhs)
{
  using size_type = typename Matrix<T, Alloc>::size_type;
  size_type ii    = 0;
  for (size_type i = 0; i < rhs.rows(); i++)
  {
    for (size_type j = 0; j < rhs.cols(); j++)
      out << rhs(ii++) << " ";
    out << std::endl;
  }
  return out;
}


template<class T, typename Alloc>
std::istream& operator>>(std::istream& is, Matrix<T, Alloc>& rhs)
{
  using size_type = typename Matrix<T, Alloc>::size_type;
  for (size_type i = 0; i < rhs.size(); i++)
  {
    is >> rhs(i++);
  }
  return is;
}
//-----------------------------------------------------------------------------
// We need to specialize CreateLeaf<T> for our class, so that operators
// know what to stick in the leaves of the expression tree.
//-----------------------------------------------------------------------------
template<class T, typename Alloc>
struct CreateLeaf<Matrix<T, Alloc>>
{
  using Leaf_t = Reference<Matrix<T, Alloc>>;
  inline static Leaf_t make(const Matrix<T, Alloc>& a) { return Leaf_t(a); }
};

//-----------------------------------------------------------------------------
// We need to write a functor that is capable of comparing the size of
// the vector with a stored value. Then, we supply LeafFunctor specializations
// for Scalar<T> and STL vector leaves.
//-----------------------------------------------------------------------------
class SizeLeaf2
{
public:
  using size_type = int;

  SizeLeaf2(size_type s, size_type p) : size_m(s), size_n(p) {}
  SizeLeaf2(const SizeLeaf2& model) : size_m(model.size_m), size_n(model.size_n) {}

  bool operator()(size_type s, size_type p) const { return ((size_m == s) && (size_n == p)); }

private:
  size_type size_m, size_n;
};

template<class T>
struct LeafFunctor<Scalar<T>, SizeLeaf2>
{
  using Type_t = bool;
  inline static bool apply(const Scalar<T>&, const SizeLeaf2&)
  {
    // Scalars always conform.
    return true;
  }
};

template<class T, typename Alloc>
struct LeafFunctor<Matrix<T, Alloc>, SizeLeaf2>
{
  using Type_t = bool;
  inline static bool apply(const Matrix<T, Alloc>& v, const SizeLeaf2& s) { return s(v.rows(), v.cols()); }
};

//-----------------------------------------------------------------------------
// EvalLeaf1 is used to evaluate expression with matrices.
// (It's already defined for Scalar values.)
//-----------------------------------------------------------------------------
//  template<class T, typename Alloc>
//  struct LeafFunctor<Matrix<T,Alloc>,EvalLeaf1>
//  {
//    using Type_t = T;
//    inline static
//    Type_t apply(const Matrix<T,Alloc>& mat, const EvalLeaf1 &f)
//    {
//      return vec[f.val1()];
//    }
//  };
//-----------------------------------------------------------------------------
// EvalLeaf2 is used to evaluate expression with matrices.
// (It's already defined for Scalar values.)
//-----------------------------------------------------------------------------
template<class T, typename Alloc>
struct LeafFunctor<Matrix<T, Alloc>, EvalLeaf2>
{
  using Type_t = T;
  inline static Type_t apply(const Matrix<T, Alloc>& mat, const EvalLeaf2& f) { return mat(f.val1(), f.val2()); }
};


///////////////////////////////////////////////////////////////////////////////
// LOOP is done by evaluate function
///////////////////////////////////////////////////////////////////////////////
template<class T, typename Alloc, class Op, class RHS>
inline void evaluate(Matrix<T, Alloc>& lhs, const Op& op, const Expression<RHS>& rhs)
{
  if (forEach(rhs, SizeLeaf2(lhs.rows(), lhs.cols()), AndCombine()))
  {
    // We get here if the vectors on the RHS are the same size as those on
    // the LHS.
    for (int i = 0; i < lhs.rows(); ++i)
    {
      for (int j = 0; j < lhs.cols(); ++j)
      {
        op(lhs(i, j), forEach(rhs, EvalLeaf2(i, j), OpCombine()));
      }
    }
  }
  else
  {
    throw std::runtime_error("Error in evaluate: LHS and RHS don't conform in OhmmsMatrix.");
  }
}
} // namespace qmcplusplus

#include "OhmmsPETE/OhmmsMatrixOperators.h"

#endif // OHMMS_PETE_MATRIX_H
