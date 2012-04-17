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
#ifndef OHMMS_PETE_ARRAY_H
#define OHMMS_PETE_ARRAY_H

/** A D-dimensional Array class based on PETE
 *
 *  Equivalent to blitz::Array<T,D>, pooma::Array<D,T>. 
 *  No operators are provided.
 *  \todo PETE
 */
#include <vector>
#include <iostream>
using namespace std;

template<class T, unsigned D>
class Array {
public:

  typedef T          Type_t;
  typedef vector<T>  Container_t;
  typedef Array<T,D> This_t;

  Array();

  Array(const Array<T,D>& rhs);// copy contructor

  // specialized for 1-Dim
  Array(unsigned n) { 
    resize(n);
  }

  // specialized for 2-Dim
  Array(unsigned n, unsigned m) {
    resize(n,m);
  }

  // specialized for 3-Dim
  Array(unsigned l, unsigned m, unsigned n) {
    resize(l,m,n);
  }

  // specialized for 4-Dim
  Array(unsigned l, unsigned m, unsigned n, unsigned o) {
    resize(l,m,n,o);
  }

  // Destructor 
  ~Array() {}

  inline unsigned size() const { return X.size();}
  int size(int i) const { return Length[i];}

  Container_t& storage()
  {
    return X;
  }

  template<typename ST1>
  void resize(ST1* newdims)
  {
    int ntot=1;
    for(int i=0; i<D; ++i) ntot *=Length[i]=newdims[i];
    if(ntot==0) return;
    X.resize(ntot);
  }

  // resize is specialized for D-dimensional
  void resize(unsigned );          
  void resize(unsigned, unsigned);  
  void resize(unsigned, unsigned, unsigned); 
  void resize(unsigned, unsigned, unsigned, unsigned);

  inline typename Container_t::iterator begin() { return X.begin();}
  inline typename Container_t::iterator end() { return X.end();}
  inline typename Container_t::const_iterator begin() const { return X.begin();}
  inline typename Container_t::const_iterator end() const { return X.end();}

  inline Type_t* data() { 
    return &(X[0]);
  }

  inline const Type_t* data() const { 
    return &(X[0]);
  }

  inline const Type_t* first_address() const { 
    return &(X[0]);
  }

  inline const Type_t* last_address() const { 
    return &(X[0])+X.size();
  }

  inline Type_t* first_address() { 
    return &(X[0]);
  }

  inline Type_t* last_address() { 
    return &(X[0])+X.size();
  }

  This_t& operator=(const T& rhs) {
    std::fill(X.begin(),X.end(),rhs);
    return *this;
  }

  // Get and Set Operations
  inline Type_t& operator()(unsigned int i) { 
    return X[i];
  }

  inline Type_t operator()(unsigned int i) const { 
    return X[i];
  }
  inline Type_t& operator() (unsigned int i, unsigned int j) { 
    return X[j+Length[1]*i];
  }
  inline Type_t operator() (unsigned int i, unsigned int j) const { 
    return X[j+Length[1]*i];
  }
  inline Type_t& operator()(unsigned int i, unsigned int j, unsigned int k) { 
    return X[k+Length[2]*(j+Length[1]*i)];
  }
  inline Type_t operator()(unsigned int i, unsigned int j, unsigned int k) const { 
    return X[k+Length[2]*(j+Length[1]*i)];
  }
  inline Type_t& operator()(unsigned int i, unsigned int j, 
                            unsigned int k, unsigned int l) { 
    return X[l+Length[3]*(k+Length[2]*(j+Length[1]*i))];
  }
  inline Type_t operator() (unsigned int i, unsigned int j, 
                            unsigned int k, unsigned int l) const { 
    return X[l+Length[3]*(k+Length[2]*(j+Length[1]*i))];
  }

private:
  int Length[D];
  Container_t X;
};


template<class T, unsigned D>
Array<T,D>::Array() {
  for(int i=0; i<D; i++) Length[i] = 0;
}

template<class T, unsigned D>
Array<T,D>::Array(const Array<T,D>& rhs) {
  // resize the matrix
  resize(rhs.X.size());
  // assign the D-dimension indices
  for(int i=0; i<D; i++) Length[i] = rhs.Length[i];
  std::copy(rhs.begin(),rhs.end(),X.begin());
}

template<class T, unsigned D>
void Array<T,D>::resize(unsigned n) { 
  Length[0] = n;
  X.resize(n,T());
}

template<class T, unsigned D>
void Array<T,D>::resize(unsigned n, unsigned m) {
  Length[0] = n;
  Length[1] = m;
  X.resize(n*m,T());
}

template<class T, unsigned D>
void Array<T,D>::resize(unsigned l, unsigned m, unsigned n) {
  Length[0] = l;
  Length[1] = m;
  Length[2] = n;
  X.resize(l*m*n,T());
}
#endif //OHMMS_PETE_ARRAY_H

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

