//////////////////////////////////////////////////////////////////
// (c) Copyright 1998-2002 by Jeongnim Kim
//
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//   Department of Physics, Ohio State University
//   Ohio Supercomputer Center
//////////////////////////////////////////////////////////////////
// -*- C++ -*-

#ifndef OHMMS_PETE_MATRIX_H
#define OHMMS_PETE_MATRIX_H

#include "PETE/PETE.h"

#include <vector>
#include <iostream>
using namespace std;

template<class T, class C = vector<T> >
class Matrix {
public:

  typedef T           Type_t;
  typedef T           value_type;
  typedef C           Container_t;
  typedef Matrix<T,C> This_t;

  Matrix():D1(0),D2(0){ } // Default Constructor initializes to zero.

  Matrix(unsigned n) { 
    resize(n,n);
    assign(*this, T());
  }

  Matrix(unsigned n, unsigned m) { 
    resize(n,m);
    assign(*this, T());
  }

  // Copy Constructor 
  Matrix(const Matrix<T,C> &rhs){
    resize(rhs.D1, rhs.D2);
    assign(*this, rhs);
  }

  // Destructor 
  ~Matrix() { }

  inline unsigned size() const { return X.size();}
  inline unsigned size(int i) const { return (i == 0)? D1: D2;}
  inline unsigned extent(int i) const { return (i == 0)? D1: D2;}
  inline unsigned rows() const { return D1;}
  inline unsigned cols() const { return D2;}

//   inline const T* begin() const { return X.begin();}
//   inline T* begin() { return X.begin();}
//   inline const T* end() const { return X.end();}
//   inline T* end()   { return X.end();}

  inline typename Container_t::iterator begin() { return X.begin();}
  inline typename Container_t::iterator end() { return X.end();}
  inline typename Container_t::const_iterator begin() const { return X.begin();}
  inline typename Container_t::const_iterator end() const { return X.end();}

  void resize(unsigned n, unsigned m); // resize to n x m
  void add(unsigned n); // you can add rows: adding columns are forbidden

  // Assignment Operators
  This_t& operator=(const Matrix<T,C> &rhs) {
    return assign(*this,rhs);
  }

  const This_t &operator=(const Matrix<T, C> &rhs) const {
    return assign(*this, rhs);
  }

  template<class RHS>
  This_t& operator=(const RHS& rhs) {
    return assign(*this,rhs);
  }

  // Get and Set Operations for assignment operators

  // returns a pointer of i-th row 
  inline Type_t* data() { 
    return &(X[0]);
  }

  // returns a pointer of i-th row 
  inline const Type_t* data() const { 
    return &(X[0]);
  }

  // returns a const pointer of i-th row 
  inline const Type_t* operator[](unsigned int i) const { 
    return &(X[0]) + i*D2;
  }

  /// returns a pointer of i-th row, g++ iterator problem
  inline Type_t* operator[](unsigned int i) { 
    return &(X[0]) + i*D2;
  }

  inline Type_t& operator()(unsigned int i) { 
    return X[i];
  }
  // returns the i-th value in D1*D2 vector
  inline Type_t operator()(unsigned int i) const { 
    return X[i];
  }

  // returns val(i,j)
  inline Type_t& operator()(unsigned int i, unsigned int j) { 
    return X[i*D2+j];
  }

  // returns val(i,j)
  inline Type_t operator()( unsigned int i, unsigned int j) const { 
    return X[i*D2+j];
  }

  //! \param sub an input array to be copied
  //! \param d1  row-dimension of the input array
  //! \param d2  column-dimension of the input array
  //! \param i0  starting row where the copying is done
  //! \param j0  starting column where the copying is done
  template<class T1>
  inline void add(const T1* sub, unsigned int d1, unsigned d2, 
                   unsigned int i0, unsigned int j0) {
    int ii=0;
    for(int i=0; i<d1; i++) {
      int kk = (i0+i)*D2 + j0;
      for(int j=0; j<d2; j++) {
        X[kk++] += sub[ii++];
      }
    }
  }

  template<class T1>
  inline void add(const T1* sub, unsigned int d1, unsigned d2, 
                  unsigned int i0, unsigned int j0, const T& phi) {
    int ii=0;
    for(int i=0; i<d1; i++) {
      int kk = (i0+i)*D2 + j0;
      for(int j=0; j<d2; j++) {
        X[kk++] += phi*sub[ii++];
      }
    }
  }
  template<class SubMat_t>
  inline void add(const SubMat_t& sub, unsigned int i0, unsigned int j0) {
    int ii=0;
    for(int i=0; i<sub.rows(); i++) {
      int kk = (i0+i)*D2 + j0;
      for(int j=0; j<sub.cols(); j++) {
        X[kk++] += sub(ii++);
      }
    }
  }
  inline void add(const This_t& sub, unsigned int i0, unsigned int j0) {
    int ii=0;
    for(int i=0; i<sub.rows(); i++) {
      int kk = (i0+i)*D2 + j0;
      for(int j=0; j<sub.cols(); j++) {
        X[kk++] += sub[ii++];
      }
    }
  }

protected:
  int D1, D2;
  Container_t X;
};

template<class T, class C>
void Matrix<T,C>::resize(unsigned n, unsigned m){
  D1 = n; D2 = m;
  X = C(n*m,T());
}

template<class T, class C>
void Matrix<T,C>::add(unsigned n){
  X.insert(X.end(), n*D2, T());
  D1 += n;  
}


// I/O
template<class T, class C>
ostream& operator<<(ostream& out, const Matrix<T,C>& rhs)
{
  int ii=0;
  for(int i=0; i<rhs.rows(); i++) {
    for(int j=0; j<rhs.cols(); j++)
      out << rhs(ii++) << " ";
    out << endl;
  }
  return out;
}


template<class T, class C>
std::istream& operator>>(std::istream& is, Matrix<T,C>& rhs)
{
  int ii=0;
  for(int i=0; i<rhs.size(); i++) {is >> rhs(i++);}
  return is;
}
//-----------------------------------------------------------------------------
// We need to specialize CreateLeaf<T> for our class, so that operators
// know what to stick in the leaves of the expression tree.
//-----------------------------------------------------------------------------
template<class T, class C>
struct CreateLeaf<Matrix<T, C> >
{
  typedef Reference<Matrix<T, C> > Leaf_t;
  inline static
  Leaf_t make(const Matrix<T, C> &a) { return Leaf_t(a); }
};

//-----------------------------------------------------------------------------
// We need to write a functor that is capable of comparing the size of
// the vector with a stored value. Then, we supply LeafFunctor specializations
// for Scalar<T> and STL vector leaves.
//-----------------------------------------------------------------------------
class SizeLeaf2
{
public:

  SizeLeaf2(int s, int p) : size_m(s), size_n(p) { }
  SizeLeaf2(const SizeLeaf2 &model) 
  : size_m(model.size_m), size_n(model.size_n) { }

  bool operator()(int s, int p) const { 
    return ((size_m == s) && (size_n ==p)); 
  }
  
private:
  
  int size_m, size_n;
  
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

template<class T, class C>
struct LeafFunctor<Matrix<T, C>, SizeLeaf2>
{
  typedef bool Type_t;
  inline static
  bool apply(const Matrix<T, C> &v, const SizeLeaf2 &s) 
  {
    return s(v.rows(), v.cols());
  }
};

//-----------------------------------------------------------------------------
// EvalLeaf1 is used to evaluate expression with matrices.
// (It's already defined for Scalar values.)
//-----------------------------------------------------------------------------
//  template<class T, class C>
//  struct LeafFunctor<Matrix<T, C>,EvalLeaf1>
//  {
//    typedef T Type_t;
//    inline static
//    Type_t apply(const Matrix<T, C>& mat, const EvalLeaf1 &f)
//    {
//      return vec[f.val1()];
//    }
//  };
//-----------------------------------------------------------------------------
// EvalLeaf2 is used to evaluate expression with matrices.
// (It's already defined for Scalar values.)
//-----------------------------------------------------------------------------
template<class T, class C>
struct LeafFunctor<Matrix<T, C>,EvalLeaf2>
{
  typedef T Type_t;
  inline static
  Type_t apply(const Matrix<T, C>& mat, const EvalLeaf2 &f)
  {
    return mat(f.val1(), f.val2());
  }
};



///////////////////////////////////////////////////////////////////////////////
// LOOP is done by evaluate function
///////////////////////////////////////////////////////////////////////////////
template<class T, class C, class Op, class RHS>
inline void evaluate(Matrix<T, C> &lhs, const Op &op, 
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
      cerr << "Error: LHS and RHS don't conform." << endl;
      exit(1);
    }
}

#include "OhmmsPETE/OhmmsMatrixOperators.h"

#endif // OHMMS_PETE_MATRIX_H

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
