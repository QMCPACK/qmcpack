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
#ifndef OHMMS_RANDOMFUNCTORS_H
#define OHMMS_RANDOMFUNCTORS_H

//////////////////////////////////////////////////////////
// generate a random position vector ([0,1), [0,1) ...)
//////////////////////////////////////////////////////////
template<class Vec, class RNG>
class RandomUniformPos {};

template<class RNG>
class RandomUniformPos<TinyVector<double,3>,RNG> {
public:
  typedef TinyVector<double,3> Return_t;
  RandomUniformPos(RNG& rg): d_engine(rg){ }
  inline Return_t operator()(){
    return Return_t(d_engine(), d_engine(), d_engine());
  }
private:
  RNG& d_engine;
};


template<class RNG>
class RandomUniformPos<TinyVector<double,2>,RNG> {
public:
  typedef TinyVector<double,2> Return_t;
  RandomUniformPos(RNG& rg): d_engine(rg){ }
  inline Return_t operator()(){
    return Return_t(d_engine(), d_engine());
  }
private:
  RNG& d_engine;
};

//////////////////////////////////////////////////////////
// generic getNormSeq functors to produce a normalized random vector
// ceterned at zero
//////////////////////////////////////////////////////////
template<class VT, class RNG>
struct NormRandomSeq { };

// specialized for TinyVector<T,D>
template<class T, unsigned D, class RNG>
struct NormRandomSeq<TinyVector<T,D>,RNG > {
  typedef TinyVector<double,D> Return_t;
  static inline Return_t get(RNG& rng) {
    Return_t res;
    for(int i=0; i<D; i++) res[i] = rng();
    T norm = 1.0/sqrt(dot(res,res));
    return norm*res;
  }
};

//
// specialized for TinyVector<double,3>
// 
template<class RNG>
struct NormRandomSeq<TinyVector<double,3>, RNG> {
  typedef TinyVector<double,3> Return_t;
  static inline Return_t get(RNG& rng) {
    double x,y,s;
    s = 2.0e0;
    while(s > 1.0e0) {
      x = 2.0e0*rng() - 1.0e0;  // [-0.5,0.5)
      y = 2.0e0*rng() - 1.0e0;  // [-0.5,0.5)
     s = x*x + y*y;
    }
    double z = 2.0e0*sqrt(1.0e0-s);
    return Return_t(z*x, z*y, 1.0e0-2.0e0*s);
  }
};


//////////////////////////////////////////////////////////
// generic RandomSeq 
// ceterned at zero [-0.5,0.5)
//////////////////////////////////////////////////////////
template<class VT, class RNG>
struct RandomSeq { };

// specialized for TinyVector<T,D>
template<class T, unsigned D, class RNG>
struct RandomSeq<TinyVector<T,D>,RNG > {
  typedef TinyVector<double,D> Return_t;
  static inline Return_t get(RNG& rng) {
    Return_t res(0);
    for(int i=0; i<D; i++) res[i] = 2*rng()-1;
    return res;
  }
};

// specialized for TinyVector<double,3>
template<class RNG>
struct RandomSeq<TinyVector<double,3>, RNG> {
  typedef TinyVector<double,3> Return_t;
  static inline Return_t get(RNG& rng) {
    return Return_t(2*rng()-1, 2*rng()-1, 2*rng()-1); 
  }
};

//  //////////////////////////////////////////////////////////
//  // Templated function to generate random Vector<T>
//  //////////////////////////////////////////////////////////

//  //
//  // Vector of normalized random avg[a] = 0
//  //
//  template<class VT, class RNG>
//  void generateNormRandomSeq(Vector<VT>& a, RNG& rng) {
//    for(int i=0; i<a.size(); i++)
//      a[i] = NormRandomSeq<VT,RNG>::get(rng);
//  }

//  //
//  // Vector of random
//  //
//  template<class VT, class RNG>
//  void generateRandomSeq(Vector<VT>& a, RNG& rng) {
//    for(int i=0; i<a.size(); i++)
//      a[i] = RandomSeq<VT,RNG>::get(rng);
//  }

//  //
//  // Specialized for double [-0.5,0.5)
//  //
//  template<class RNG>
//  void generateRandomSeq(Vector<double>& a, RNG& rng) {
//    for(int i=0; i<a.size(); i++) a[i] = 2.0*rng()-1;
//  }


// generic templated class to generate a random vector
template<class VT, class RNG>
struct RandomVector { };


// generic templated class to generate a random sequence
template<class RA, class RNG>
struct RandomSequence { };

// specialized for vector<T> to generate a random sequence of type T
template<class T, class RNG>
struct RandomSequence<vector<T>, RNG> { 

  static void apply(vector<T>& v, RandomVector<T,RNG>& rnd) {
    typename vector<T>::iterator it= v.begin();
    while(it != v.end()) {
      (*it) = rnd(); it++;
    }
  }
};

// specialized for vector<double> to generate a random sequence of type T
template<class RNG>
struct RandomSequence<vector<double>, RNG> {

  static void apply(vector<double>& s, RNG& rnd) {
    typename vector<double>::iterator it= s.begin();
    while(it != s.end()) {
      (*it) = rnd(); it++;
    }
  }
};


#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
