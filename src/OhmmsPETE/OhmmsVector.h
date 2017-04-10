//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//		      Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign   
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign 
//////////////////////////////////////////////////////////////////////////////////////


#ifndef OHMMS_PETE_VECTOR_H
#define OHMMS_PETE_VECTOR_H

/** A one-dimensional vector class based on PETE
 *
 *  Closely related to PETE STL vector example.
 *  Equivalent to blitz::Array<T,1>, pooma::Array<1,T>.
 *  class C is a container class. Default is std::vector<T>
 *  \todo Implement openMP compatible container class or evaluate function.
 *  \todo Implement get/put member functions for MPI-like parallelism
 */


#include "PETE/PETE.h"

#include <cstdlib>
#include <vector>
#include <iostream>
//
namespace qmcplusplus
{
template<class T, class C = std::vector<T> >
struct Vector
{

  typedef T           Type_t;
  typedef C           Container_t;
  typedef Vector<T,C> This_t;
  typedef typename Container_t::iterator iterator;
  typedef typename Container_t::const_iterator const_iterator;

  Vector() { } // Default Constructor initializes to zero.

  Vector(unsigned n)
  {
    resize(n);
    assign(*this, T());
  }

  // Copy Constructor
  Vector(const Vector<T,C> &rhs)
  {
    resize(rhs.size());
    assign(*this,rhs);
  }

  // Destructor
  ~Vector() { }

  inline unsigned size() const
  {
    return X.size();
  }
  void resize(unsigned n); // resize to n
  void create(unsigned n); // add n elements

  // Assignment Operators
  This_t& operator=(const Vector<T,C> &rhs)
  {
    return assign(*this,rhs);
  }

  const This_t &operator=(const Vector<T, C> &rhs) const
  {
    return assign(*this, rhs);
  }

  template<class RHS>
  This_t& operator=(const RHS& rhs)
  {
    return assign(*this,rhs);
  }

  inline Type_t* data()
  {
    return &(X[0]);
  }
  inline const Type_t* data() const
  {
    return &(X[0]);
  }

  inline Type_t* first_address()
  {
    return &(X[0]);
  }
  inline const Type_t* first_address() const
  {
    return &(X[0]);
  }

  inline Type_t* last_address()
  {
    return &(X[0])+size();
  }
  inline const Type_t* last_address() const
  {
    return &(X[0])+size();
  }

  inline iterator begin()
  {
    return X.begin();
  }
  inline const_iterator begin() const
  {
    return X.begin();
  }

  inline iterator end()
  {
    return X.end();
  }
  inline const_iterator end() const
  {
    return X.end();
  }

  //inline Type_t* begin() { return X.begin();}
  //inline const Type_t* begin() const { return X.begin();}

  //inline Type_t* end() { return X.end();}
  //inline const Type_t* end() const { return X.end();}

  // Get and Set Operations
  inline Type_t& operator[](unsigned int i)
  {
    return X[i];
  }

  inline Type_t operator[](unsigned int i) const
  {
    return X[i];
  }

  inline Type_t& operator()(unsigned int i)
  {
    return X[i];
  }

  inline Type_t operator()( unsigned int i) const
  {
    return X[i];
  }

  //----------------------------------------------------------------------
  // parallel communication

  //Message& putMessage(Message& m) const {
  //  m.setCopy(true);
  //  ::putMessage(m, X, X + D);
  //    return m;
  //}

  //Message& getMessage(Message& m) {
  //  ::getMessage(m, X, X + D);
  //  return m;
  //}


private:
  Container_t X;
};

template<class T, class C>
void Vector<T,C>::resize(unsigned n)
{
  X = C(n,T());
}

template<class T, class C>
void Vector<T,C>::create(unsigned n)
{
  X.insert(X.end(), n, T());
}
}

#include "OhmmsPETE/OhmmsVectorOperators.h"

namespace qmcplusplus
{
//-----------------------------------------------------------------------------
// We need to specialize CreateLeaf<T> for our class, so that operators
// know what to stick in the leaves of the expression tree.
//-----------------------------------------------------------------------------
template<class T, class C>
struct CreateLeaf<Vector<T, C> >
{
  typedef Reference<Vector<T, C> > Leaf_t;
  inline static
  Leaf_t make(const Vector<T, C> &a)
  {
    return Leaf_t(a);
  }
};

//-----------------------------------------------------------------------------
// We need to write a functor that is capable of comparing the size of
// the vector with a stored value. Then, we supply LeafFunctor specializations
// for Scalar<T> and STL vector leaves.
//-----------------------------------------------------------------------------
class SizeLeaf
{
public:

  SizeLeaf(int s) : size_m(s) { }
  SizeLeaf(const SizeLeaf &model) : size_m(model.size_m) { }
  bool operator()(int s) const
  {
    return size_m == s;
  }

private:

  int size_m;

};

template<class T>
struct LeafFunctor<Scalar<T>, SizeLeaf>
{
  typedef bool Type_t;
  inline static
  bool apply(const Scalar<T> &, const SizeLeaf &)
  {
    // Scalars always conform.
    return true;
  }
};

template<class T, class C>
struct LeafFunctor<Vector<T, C>, SizeLeaf>
{
  typedef bool Type_t;
  inline static
  bool apply(const Vector<T, C> &v, const SizeLeaf &s)
  {
    return s(v.size());
  }
};

//-----------------------------------------------------------------------------
// EvalLeaf1 is used to evaluate expression with vectors.
// (It's already defined for Scalar values.)
//-----------------------------------------------------------------------------

template<class T, class C>
struct LeafFunctor<Vector<T, C>,EvalLeaf1>
{
  typedef T Type_t;
  inline static
  Type_t apply(const Vector<T, C>& vec,const EvalLeaf1 &f)
  {
    return vec[f.val1()];
  }
};

//////////////////////////////////////////////////////////////////////////////////
// LOOP is done by evaluate function
//////////////////////////////////////////////////////////////////////////////////
template<class T, class C, class Op, class RHS>
inline void evaluate(Vector<T, C> &lhs, const Op &op, const Expression<RHS> &rhs)
{
  if (forEach(rhs, SizeLeaf(lhs.size()), AndCombine()))
  {
    // We get here if the vectors on the RHS are the same size as those on
    // the LHS.
    for (int i = 0; i < lhs.size(); ++i)
    {
      // The actual assignment operation is performed here.
      // PETE operator tags all define operator() to perform the operation.
      // (In this case op performs an assignment.) forEach is used
      // to compute the rhs value.  EvalLeaf1 gets the
      // values at each node using random access, and the tag
      // OpCombine tells forEach to use the operator tags in the expression
      // to combine values together.
      op(lhs[i], forEach(rhs, EvalLeaf1(i), OpCombine()));
    }
  }
  else
  {
    std::cerr << "Error: LHS and RHS don't conform in OhmmsVector." << std::endl;
    abort();
  }
}
// I/O
template<class T, class C>
std::ostream& operator<<(std::ostream& out, const Vector<T,C>& rhs)
{
  for(int i=0; i<rhs.size(); i++)
    out << rhs[i] << std::endl;
  return out;
}

template<class T, class C>
std::istream& operator>>(std::istream& is, Vector<T,C>& rhs)
{
  //printTinyVector<TinyVector<T,D> >::print(out,rhs);
  for(int i=0; i<rhs.size(); i++)
    is >> rhs[i];
  return is;
}

}
#endif // VEKTOR_H

