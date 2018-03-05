//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//		      Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//  		      Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign   
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign 
//////////////////////////////////////////////////////////////////////////////////////

#ifndef OHMMS_PETE_MATRIX_H
#define OHMMS_PETE_MATRIX_H

#include "PETE/PETE.h"
#include <cstdlib>
#include "OhmmsPETE/OhmmsVector.h"
#include <iostream>

namespace qmcplusplus
{

template<class T, typename Alloc=std::allocator<T> >
class Matrix
{
public:

  typedef T            Type_t;
  typedef T            value_type;
  typedef T*           pointer;
  typedef const T*     const_pointer;
  typedef Vector<T,Alloc>  Container_t;
  typedef typename Container_t::size_type size_type;
  typedef typename Container_t::iterator iterator;
  typedef Matrix<T,Alloc>  This_t;

  Matrix():D1(0),D2(0),TotSize(0) { } // Default Constructor initializes to zero.

  Matrix(size_type n)
  {
    resize(n,n);
    //assign(*this, T());
  }

  Matrix(size_type n, size_type m)
  {
    resize(n,m);
    //assign(*this, T());
  }

  // Copy Constructor
  Matrix(const Matrix<T,Alloc> &rhs)
  {
    copy(rhs);
  }

  // Destructor
  ~Matrix() { }

  inline size_type size() const
  {
    return TotSize;
  }
  inline size_type rows() const
  {
    return D1;
  }
  inline size_type cols() const
  {
    return D2;
  }
  inline size_type size1() const
  {
    return D1;
  }
  inline size_type size2() const
  {
    return D2;
  }
  inline size_type size(int i) const
  {
    return (i == 0)? D1: D2;
  }
  inline size_type extent(int i) const
  {
    return (i == 0)? D1: D2;
  }

//   inline const T* begin() const { return X.begin();}
//   inline T* begin() { return X.begin();}
//   inline const T* end() const { return X.end();}
//   inline T* end()   { return X.end();}

  inline typename Container_t::iterator begin()
  {
    return X.begin();
  }
  inline typename Container_t::iterator end()
  {
    return X.end();
  }
  inline typename Container_t::const_iterator begin() const
  {
    return X.begin();
  }
  inline typename Container_t::const_iterator end() const
  {
    return X.end();
  }

  inline typename Container_t::iterator begin(int i)
  {
    return X.begin()+i*D2;
  }
  inline typename Container_t::const_iterator begin(int i) const
  {
    return X.begin()+i*D2;
  }

  inline void resize(size_type n, size_type m)
  {
    D1 = n;
    D2 = m;
    TotSize=n*m;
    X.resize(n*m);
  }

  // free the matrix storage
  inline void free()
  {
    X.free();
  }

  // Attach to pre-allocated memory
  inline void attachReference(T* ref)
  {
    X.attachReference(ref, TotSize);
  }

  inline void add(size_type n)   // you can add rows: adding columns are forbidden
  {
    X.insert(X.end(), n*D2, T());
    D1 += n;
  }

  inline void copy(const Matrix<T,Alloc>& rhs)
  {
    resize(rhs.D1, rhs.D2);
    assign(*this, rhs);
  }

  // Assignment Operators
  inline This_t& operator=(const Matrix<T,Alloc> &rhs)
  {
    resize(rhs.D1,rhs.D2);
    return assign(*this,rhs);
  }

  inline const This_t &operator=(const Matrix<T,Alloc> &rhs) const
  {
    return assign(*this, rhs);
  }

  template<class RHS>
  This_t& operator=(const RHS& rhs)
  {
    return assign(*this,rhs);
  }

  // Get and Set Operations for assignment operators
  // returns a pointer of i-th row
  inline pointer data()
  {
    return X.data();
  }

  // returns a pointer of i-th row
  inline const_pointer data() const
  {
    return X.data();
  }

  // returns a const pointer of i-th row
  inline const Type_t* data(size_type i) const
  {
    return X.data() + i*D2;
  }

  /// returns a pointer of i-th row, g++ iterator problem
  inline Type_t* data(size_type i)
  {
    return X.data() + i*D2;
  }

  inline pointer first_address()
  {
    return X.data();
  }

  // returns a pointer of i-th row
  inline const_pointer first_address() const
  {
    return X.data();
  }

  inline pointer last_address()
  {
    return X.data()+TotSize;
  }

  // returns a pointer of i-th row
  inline const Type_t* last_address() const
  {
    return X.data()+TotSize;
  }


  // returns a const pointer of i-th row
  inline const Type_t* operator[](size_type i) const
  {
    return X.data() + i*D2;
  }

  /// returns a pointer of i-th row, g++ iterator problem
  inline Type_t* operator[](size_type i)
  {
    return X.data() + i*D2;
  }

  inline Type_t& operator()(size_type i)
  {
    return X[i];
  }
  // returns the i-th value in D1*D2 vector
  inline Type_t operator()(size_type i) const
  {
    return X[i];
  }

  // returns val(i,j)
  inline Type_t& operator()(size_type i, size_type j)
  {
    return X[i*D2+j];
  }

  // returns val(i,j)
  inline Type_t operator()( size_type i, size_type j) const
  {
    return X[i*D2+j];
  }

  inline void swap_rows (int r1, int r2)
  {
    for (int col=0; col<D2; col++)
    {
      Type_t tmp = (*this)(r1,col);
      (*this)(r1,col) = (*this)(r2,col);
      (*this)(r2,col) = tmp;
    }
  }

  inline void swap_cols (int c1, int c2)
  {
    for (int row=0; row<D1; row++)
    {
      Type_t tmp = (*this)(row, c1);
      (*this)(row, c1) = (*this)(row, c2);
      (*this)(row, c2) = tmp;
    }
  }


  template<class IT>
  inline void replaceRow(IT first, size_type i)
  {
    std::copy(first,first+D2,X.begin()+i*D2);
  }

  template<class IT>
  inline void replaceColumn(IT first,size_type j)
  {
    typename Container_t::iterator ii(X.begin()+j);
    for(int i=0; i<D1; i++, ii+=D2)
      *ii=*first++;
  }

  template<class IT>
  inline void add2Column(IT first,size_type j)
  {
    typename Container_t::iterator ii(X.begin()+j);
    for(int i=0; i<D1; i++, ii+=D2)
      *ii+=*first++;
  }

  /**
   * \param sub an input array to be copied
   * \param d1  row-dimension of the input array
   * \param d2  column-dimension of the input array
   * \param i0  starting row where the copying is done
   * \param j0  starting column where the copying is done
   */
  template<class T1>
  inline void add(const T1* sub, size_type d1, size_type d2, size_type i0, size_type j0)
  {
    int ii=0;
    for(int i=0; i<d1; i++)
    {
      int kk = (i0+i)*D2 + j0;
      for(int j=0; j<d2; j++)
      {
        X[kk++] += sub[ii++];
      }
    }
  }

  template<class T1>
  inline void add(const T1* sub, size_type d1, size_type d2, size_type i0, size_type j0, const T& phi)
  {
    size_type ii=0;
    for(size_type i=0; i<d1; i++)
    {
      int kk = (i0+i)*D2 + j0;
      for(size_type j=0; j<d2; j++)
      {
        X[kk++] += phi*sub[ii++];
      }
    }
  }
  template<class SubMat_t>
  inline void add(const SubMat_t& sub, unsigned int i0, unsigned int j0)
  {
    size_type ii=0;
    for(size_type i=0; i<sub.rows(); i++)
    {
      int kk = (i0+i)*D2 + j0;
      for(size_type j=0; j<sub.cols(); j++)
      {
        X[kk++] += sub(ii++);
      }
    }
  }
  inline void add(const This_t& sub, unsigned int i0, unsigned int j0)
  {
    size_type ii=0;
    for(size_type i=0; i<sub.rows(); i++)
    {
      int kk = (i0+i)*D2 + j0;
      for(size_type j=0; j<sub.cols(); j++)
      {
        X[kk++] += sub[ii++];
      }
    }
  }

  template<class Msg>
  inline Msg& putMessage(Msg& m)
  {
    m.Pack(X.data(),D1*D2);
    return m;
  }

  template<class Msg>
  inline Msg& getMessage(Msg& m)
  {
    m.Unpack(X.data(),D1*D2);
    return m;
  }

protected:
  size_type D1, D2;
  size_type TotSize;
  Container_t X;
};

// I/O
template<class T, typename Alloc>
std::ostream& operator<<(std::ostream& out, const Matrix<T,Alloc>& rhs)
{
  typedef typename Matrix<T,Alloc>::size_type size_type;
  size_type ii=0;
  for(size_type i=0; i<rhs.rows(); i++)
  {
    for(size_type j=0; j<rhs.cols(); j++)
      out << rhs(ii++) << " ";
    out << std::endl;
  }
  return out;
}


template<class T, typename Alloc>
std::istream& operator>>(std::istream& is, Matrix<T,Alloc>& rhs)
{
  typedef typename Matrix<T,Alloc>::size_type size_type;
  size_type ii=0;
  for(size_type i=0; i<rhs.size(); i++)
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
struct CreateLeaf<Matrix<T,Alloc> >
{
  typedef Reference<Matrix<T,Alloc> > Leaf_t;
  inline static
  Leaf_t make(const Matrix<T,Alloc> &a)
  {
    return Leaf_t(a);
  }
};

//-----------------------------------------------------------------------------
// We need to write a functor that is capable of comparing the size of
// the vector with a stored value. Then, we supply LeafFunctor specializations
// for Scalar<T> and STL vector leaves.
//-----------------------------------------------------------------------------
class SizeLeaf2
{
public:

  typedef int size_type;

  SizeLeaf2(size_type s, size_type p) : size_m(s), size_n(p) { }
  SizeLeaf2(const SizeLeaf2 &model)
    : size_m(model.size_m), size_n(model.size_n) { }

  bool operator()(size_type s, size_type p) const
  {
    return ((size_m == s) && (size_n ==p));
  }

private:

  size_type size_m, size_n;
};

template<class T>
struct LeafFunctor<Scalar<T>, SizeLeaf2>
{
  typedef bool Type_t;
  inline static
  bool apply(const Scalar<T> &, const SizeLeaf2 &)
  {
    // Scalars always conform.
    return true;
  }
};

template<class T, typename Alloc>
struct LeafFunctor<Matrix<T,Alloc>, SizeLeaf2>
{
  typedef bool Type_t;
  inline static
  bool apply(const Matrix<T,Alloc> &v, const SizeLeaf2 &s)
  {
    return s(v.rows(), v.cols());
  }
};

//-----------------------------------------------------------------------------
// EvalLeaf1 is used to evaluate expression with matrices.
// (It's already defined for Scalar values.)
//-----------------------------------------------------------------------------
//  template<class T, typename Alloc>
//  struct LeafFunctor<Matrix<T,Alloc>,EvalLeaf1>
//  {
//    typedef T Type_t;
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
struct LeafFunctor<Matrix<T,Alloc>,EvalLeaf2>
{
  typedef T Type_t;
  inline static
  Type_t apply(const Matrix<T,Alloc>& mat, const EvalLeaf2 &f)
  {
    return mat(f.val1(), f.val2());
  }
};



///////////////////////////////////////////////////////////////////////////////
// LOOP is done by evaluate function
///////////////////////////////////////////////////////////////////////////////
template<class T, typename Alloc, class Op, class RHS>
inline void evaluate(Matrix<T,Alloc> &lhs, const Op &op,
                     const Expression<RHS> &rhs)
{
  if (forEach(rhs, SizeLeaf2(lhs.rows(), lhs.cols()), AndCombine()))
  {
    // We get here if the vectors on the RHS are the same size as those on
    // the LHS.
    int ii=0;
    for(int i=0; i<lhs.rows(); ++i)
    {
      for (int j = 0; j < lhs.cols(); ++j)
      {
        op(lhs(ii++), forEach(rhs, EvalLeaf2(i,j), OpCombine()));
      }
    }
  }
  else
  {
    std::cerr << "Error: LHS and RHS don't conform in OhmmsMatrix." << std::endl;
    abort();
  }
}
}

#include "OhmmsPETE/OhmmsMatrixOperators.h"

#endif // OHMMS_PETE_MATRIX_H

