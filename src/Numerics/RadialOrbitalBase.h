//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



/**@file RadialOrbitalBase.h
 * @brief declare/define the base class for RadialOrbital
 */
#ifndef QMCPLUSPLUS_RADIALORBITAL_BASE_H
#define QMCPLUSPLUS_RADIALORBITAL_BASE_H

/** base class for RadialOrbital to facilitate maping between a group of radial functors to a numerical functor
 */
template<class T>
struct RadialOrbitalBase
{
  inline RadialOrbitalBase() {}
  virtual ~RadialOrbitalBase() {}
  virtual T f(T r) const = 0;
  virtual T df(T r) const = 0;
};

/** composite class that contains a number of radial functions that belong to a group.
 */
template<class T>
struct RadialOrbitalSet: public RadialOrbitalBase<T>
{

  std::vector<RadialOrbitalBase<T>*> InFunc;

  ~RadialOrbitalSet()
  {
    typename std::vector<RadialOrbitalBase<T>*>::iterator it(InFunc.begin());
    typename std::vector<RadialOrbitalBase<T>*>::iterator it_end(InFunc.end());
    while(it != it_end)
    {
      delete *it;
      ++it;
    }
  }

  inline
  void addRadialOrbital(RadialOrbitalBase<T>* arad)
  {
    InFunc.push_back(arad);
  }

  inline T f(T r) const
  {
    typename std::vector<RadialOrbitalBase<T>*>::const_iterator it(InFunc.begin());
    typename std::vector<RadialOrbitalBase<T>*>::const_iterator it_end(InFunc.end());
    T res(0.0);
    while(it != it_end)
    {
      res += (*it)->f(r);
      ++it;
    }
    return res;
  }

  inline T df(T r) const
  {
    typename std::vector<RadialOrbitalBase<T>*>::const_iterator it(InFunc.begin());
    typename std::vector<RadialOrbitalBase<T>*>::const_iterator it_end(InFunc.end());
    T res(0.0);
    while(it != it_end)
    {
      res += (*it)->df(r);
      ++it;
    }
    return res;
  }

};
#endif
