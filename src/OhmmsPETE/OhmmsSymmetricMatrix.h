//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign   
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign 
//////////////////////////////////////////////////////////////////////////////////////


#ifndef OHMMS_PETE_SYMMETRIC_MATRIX_H
#define OHMMS_PETE_SYMMETRIC_MATRIX_H

#include "PETE/PETE.h"

#include <vector>
#include <iostream>

template<class T, class C = std::vector<T> >

class SymmetricMatrix
{
public:

  typedef T           Type_t;
  typedef C           Container_t;
  typedef SymmetricMatrix<T,C> This_t;
  typedef typename Container_t::iterator iterator;
  typedef typename Container_t::const_iterator const_iterator;

  SymmetricMatrix():D1(0),D2(0) { } // Default Constructor initializes to zero.

  SymmetricMatrix(unsigned n)
  {
    resize(n,n);
    //assign(*this, T());
  }

  SymmetricMatrix(unsigned n, unsigned m)
  {
    resize(n,m);
    //assign(*this, T());
  }

  // Copy Constructor
  SymmetricMatrix(const SymmetricMatrix<T,C> &rhs)
  {
    resize(rhs.D1, rhs.D2);
    //assign(*this, rhs);
  }

  // Destructor
  ~SymmetricMatrix() { }

  inline unsigned size() const
  {
    return X.size();
  }
  inline unsigned size(int i) const
  {
    return (i == 0)? D1: D2;
  }
  inline unsigned nrows() const
  {
    return D1;
  }
  inline unsigned ncols() const
  {
    return D2;
  }

  void clear()
  {
    std::fill(X.begin(), X.end(),T());
  }

  inline const_iterator begin() const
  {
    return X.begin();
  }
  inline iterator begin()
  {
    return X.begin();
  }
  inline const_iterator end() const
  {
    return X.end();
  }

  inline const T* data() const
  {
    return &(X[0]);
  }
  inline T* data()
  {
    return &(X[0]);
  }

  void resize(unsigned n, unsigned m); // resize to n x m

  // Assignment Operators
  This_t& operator=(const SymmetricMatrix<T,C> &rhs)
  {
    X = rhs.X;
    //return assign(*this,rhs);
  }

  const This_t &operator=(const SymmetricMatrix<T, C> &rhs) const
  {
    X = rhs.X;
    return *this;
  }

  // Get and Set Operations for assignment operators
  // returns a pointer of i-th row
  inline Type_t* operator[](unsigned int i)
  {
    //return X.begin() + i*D2;
    return &(X[i*D2]);
  }

  // returns a const pointer of i-th row
  inline const Type_t* operator[](unsigned int i) const
  {
    //return X.begin() + i*D2;
    return &(X[(2*D2-i-1)*i/2]);
  }

  // returns the i-th value to be assigned in D1*D2 vector
  inline Type_t& operator()(unsigned int i)
  {
    return X[i];
  }
  // returns the i-th value in D1*D2 vector
  inline Type_t operator()(unsigned int i) const
  {
    return X[i];
  }

  // returns val(i,j)
  inline Type_t& operator()(unsigned int i, unsigned int j)
  {
    return X[(2*D2-i-1)*i/2+j];
  }

  // returns val(i,j)
  inline Type_t operator()( unsigned int i, unsigned int j) const
  {
    return X[(2*D2-i-1)*i/2+j];
  }


protected:
  int D1, D2;
  Container_t X;
};

template<class T, class C>
void SymmetricMatrix<T,C>::resize(unsigned n, unsigned m)
{
  D1 = n;
  D2 = m;
  X = C(D2*(D2+1)/2);
}


//  //-----------------------------------------------------------------------------
//  // We need to specialize CreateLeaf<T> for our class, so that operators
//  // know what to stick in the leaves of the expression tree.
//  //-----------------------------------------------------------------------------
//  template<class T, class C>
//  struct CreateLeaf<SymmetricMatrix<T, C> >
//  {
//    typedef Reference<SymmetricMatrix<T, C> > Leaf_t;
//    inline static
//    Leaf_t make(const SymmetricMatrix<T, C> &a) { return Leaf_t(a); }
//  };

//  //-----------------------------------------------------------------------------
//  // We need to write a functor that is capable of comparing the size of
//  // the vector with a stored value. Then, we supply LeafFunctor specializations
//  // for Scalar<T> and STL vector leaves.
//  //-----------------------------------------------------------------------------
//  class SizeLeaf2
//  {
//  public:

//    SizeLeaf2(int s, int p) : size_m(s), size_n(p) { }
//    SizeLeaf2(const SizeLeaf2 &model)
//    : size_m(model.size_m), size_n(model.size_n) { }

//    bool operator()(int s, int p) const {
//      return ((size_m == s) && (size_n ==p));
//    }

//  private:

//    int size_m, size_n;

//  };

//  template<class T>
//  struct LeafFunctor<Scalar<T>, SizeLeaf2>
//  {
//    typedef bool Type_t;
//    inline static
//    bool apply(const Scalar<T> &, const SizeLeaf2 &)
//    {
//      // Scalars always conform.
//      return true;
//    }
//  };

//  template<class T, class C>
//  struct LeafFunctor<SymmetricMatrix<T, C>, SizeLeaf2>
//  {
//    typedef bool Type_t;
//    inline static
//    bool apply(const SymmetricMatrix<T, C> &v, const SizeLeaf2 &s)
//    {
//      return s(v.nrows(), v.ncols());
//    }
//  };

//-----------------------------------------------------------------------------
// EvalLeaf1 is used to evaluate expression with matrices.
// (It's already defined for Scalar values.)
//-----------------------------------------------------------------------------
//  template<class T, class C>
//  struct LeafFunctor<SymmetricMatrix<T, C>,EvalLeaf1>
//  {
//    typedef T Type_t;
//    inline static
//    Type_t apply(const SymmetricMatrix<T, C>& mat, const EvalLeaf1 &f)
//    {
//      return vec[f.val1()];
//    }
//  };
//-----------------------------------------------------------------------------
// EvalLeaf2 is used to evaluate expression with matrices.
// (It's already defined for Scalar values.)
//-----------------------------------------------------------------------------
//  template<class T, class C>
//  struct LeafFunctor<SymmetricMatrix<T, C>,EvalLeaf2>
//  {
//    typedef T Type_t;
//    inline static
//    Type_t apply(const SymmetricMatrix<T, C>& mat, const EvalLeaf2 &f)
//    {
//      return mat(f.val1(), f.val2());
//    }
//  };



///////////////////////////////////////////////////////////////////////////////
// LOOP is done by evaluate function
///////////////////////////////////////////////////////////////////////////////
template<class T, class C, class Op, class RHS>
inline void evaluate(SymmetricMatrix<T, C> &lhs, const Op &op,
                     const Expression<RHS> &rhs)
{
  if (forEach(rhs, SizeLeaf2(lhs.nrows(), lhs.ncols()), AndCombine()))
  {
    // We get here if the vectors on the RHS are the same size as those on
    // the LHS.
    int ii=0;
    for(int i=0; i<lhs.nrows(); ++i)
    {
      for (int j = 0; j < lhs.ncols(); ++j)
      {
        op(lhs(ii++), forEach(rhs, EvalLeaf2(i,j), OpCombine()));
      }
    }
  }
  else
  {
    std::cerr << "Error: LHS and RHS don't conform in OhmmsSymmetricMatrix." << std::endl;
    exit(1);
  }
}

// I/O
template<class T, class C>
ostream& operator<<(std::ostream& out, const SymmetricMatrix<T,C>& rhs)
{
  int ii=0;
  for(int i=0; i<rhs.nrows(); i++)
  {
    for(int j=0; j<i; j++)
      out << rhs(j,i) << " ";
    for(int j=i; j<rhs.ncols(); j++)
      out << rhs(i,j) << " ";
    out << std::endl;
  }
  return out;
}

#endif // OHMMS_PETE_MATRIX_H

