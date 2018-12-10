//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign   
//		      Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign 
//////////////////////////////////////////////////////////////////////////////////////


#ifndef OHMMS_PETE_ARRAY_H
#define OHMMS_PETE_ARRAY_H

/** A D-dimensional Array class based on PETE
 *
 *  Equivalent to blitz::Array<T,D>, pooma::Array<D,T>.
 *  No operators are provided.
 *  \todo use aligned_vector<T>
 */
#include <vector>
#include <iostream>

template<class T, unsigned D>
struct Array
{

  typedef T          Type_t;
  typedef std::vector<T>  Container_t; 
  typedef Array<T,D> This_t;

  size_t Length[D];
  Container_t X;

  //default constructor
  Array()
  {
    for(int i=0; i<D; i++)
      Length[i] = 0;
  }

  //copy constructor
  Array(const Array& rhs)
  {
    resize(rhs);
    std::copy(rhs.begin(),rhs.end(),X.begin());
  }

  // specialized for 1-Dim
  Array(size_t n)
  {
    resize(n);
  }

  // specialized for 2-Dim
  Array(size_t n, size_t m)
  {
    resize(n,m);
  }

  // specialized for 3-Dim
  Array(size_t l, size_t m, size_t n)
  {
    resize(l,m,n);
  }

  // specialized for 4-Dim
  Array(size_t l, size_t m, size_t n, size_t o)
  {
    resize(l,m,n,o);
  }

  // specialized for something like TinyVector
  template<typename ST1>
  Array(ST1* dims)
  {
    resize(dims);
  }

  // do nothing Destructor
  ~Array() {}

  inline unsigned dim() const
  {
    return D;
  }
  inline size_t* shape() const
  {
    return &Length[0];
  }
  inline size_t size() const
  {
    return X.size();
  }
  inline size_t size(int i) const
  {
    return Length[i];
  }

  Container_t& storage()
  {
    return X;
  }

  template<typename TT>
    void resize(const Array<TT,D>& rhs)
    {
      X.resize(rhs.size());
      for(int i=0; i<D; i++)
        Length[i] = rhs.Length[i];
    }


  template<typename ST1>
  void resize(ST1* newdims)
  {
    int ntot=1;
    for(int i=0; i<D; ++i)
      ntot *=Length[i]=newdims[i];
    if(ntot==0)
      return;
    X.resize(ntot);
  }

  // resize is specialized for D-dimensional
  void resize(size_t );
  void resize(size_t, size_t);
  void resize(size_t, size_t, size_t);
  void resize(size_t, size_t, size_t, size_t);

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

  inline Type_t* data()
  {
    return &(X[0]);
  }

  inline const Type_t* data() const
  {
    return &(X[0]);
  }

  inline const Type_t* first_address() const
  {
    return &(X[0]);
  }

  inline const Type_t* last_address() const
  {
    return &(X[0])+X.size();
  }

  inline Type_t* first_address()
  {
    return &(X[0]);
  }

  inline Type_t* last_address()
  {
    return &(X[0])+X.size();
  }

  This_t& operator=(const T& rhs)
  {
    std::fill(X.begin(),X.end(),rhs);
    return *this;
  }

  This_t& operator=(const Array& rhs)
  {
    if(&rhs != this)
    {
      resize(rhs);
      std::copy(rhs.begin(),rhs.end(),X.begin());
    }
    return *this;
  }

  template<typename TT>
  This_t& operator=(const Array<TT, D>& rhs)
  {
    resize(rhs);
    std::copy(rhs.begin(),rhs.end(),X.begin());
    return *this;
  }

  // Get and Set Operations
  inline Type_t& operator()(size_t i)
  {
    return X[i];
  }

  inline Type_t operator()(size_t i) const
  {
    return X[i];
  }
  inline Type_t& operator() (size_t i, size_t j)
  {
    return X[j+Length[1]*i];
  }
  inline Type_t operator() (size_t i, size_t j) const
  {
    return X[j+Length[1]*i];
  }
  inline Type_t& operator()(size_t i, size_t j, size_t k)
  {
    return X[k+Length[2]*(j+Length[1]*i)];
  }
  inline Type_t operator()(size_t i, size_t j, size_t k) const
  {
    return X[k+Length[2]*(j+Length[1]*i)];
  }
  inline Type_t& operator()(size_t i, size_t j,
                            size_t k, size_t l)
  {
    return X[l+Length[3]*(k+Length[2]*(j+Length[1]*i))];
  }
  inline Type_t operator() (size_t i, size_t j,
                            size_t k, size_t l) const
  {
    return X[l+Length[3]*(k+Length[2]*(j+Length[1]*i))];
  }

  inline Type_t sum()
  {
    Type_t s=0;
    for(int i=0; i<X.size(); ++i)
      s+=X[i];
    return s;
  }

};


//need to assert
template<class T, unsigned D>
void Array<T,D>::resize(size_t n)
{
  Length[0] = n;
  X.resize(n,T());
}

template<class T, unsigned D>
void Array<T,D>::resize(size_t n, size_t m)
{
  Length[0] = n;
  Length[1] = m;
  X.resize(n*m,T());
}

template<class T, unsigned D>
void Array<T,D>::resize(size_t l, size_t m, size_t n)
{
  Length[0] = l;
  Length[1] = m;
  Length[2] = n;
  X.resize(l*m*n);//,T());
}

template<class T, unsigned D>
void Array<T,D>::resize(size_t l, size_t m, size_t n, size_t o)
{
  Length[0] = l;
  Length[1] = m;
  Length[2] = n;
  Length[3] = o;
  X.resize(l*m*n*o);//,T());
}
#endif //OHMMS_PETE_ARRAY_H


