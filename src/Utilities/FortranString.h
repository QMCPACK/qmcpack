//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



#ifndef OHMMS_FORTRANSTRING_CONVERSION_H
#define OHMMS_FORTRANSTRING_CONVERSION_H

#include <stdio.h>
#include <vector>
#include <string>
#include <iostream>


/*!\class template<int N> struct FortranString
 *\brief A wrapper class to pass a std::string to fortran codes
 */
template<int N>
struct FortranString
{

  char m_data[N];

  FortranString()
  {
    for(int i=0; i<N; i++)
      m_data[i] = ' ';
  }

  FortranString(const char* c)
  {
    for(int i=0; i<N; i++)
      m_data[i] = ' ';
    set(c);
  }


  void set(const char* c)
  {
    sprintf(m_data,"%s",c);
  }

  //!< returns the pointer of the first element (same interface as std::string::c_str())
  const char* c_str() const
  {
    return  m_data;
  }

};


/*!\class template<int N> struct FortranStringArray
 *\brief A wrapper class to pass an array of strings to fortran codes
 * Simply remove "\n" so that fortran character manipulations work
 */
template<int N>
struct FortranStringArray
{

  std::vector<char> m_data;

  FortranStringArray() {}

  explicit FortranStringArray(int num)
  {
    resize(num);
  }

  ~FortranStringArray() { }

  int size() const
  {
    return m_data.size()/N;
  }

  void resize(int num)
  {
    m_data.resize(num*N,' ');
  }

  void set(const std::vector<std::string>& a)
  {
    resize(a.size());
    for(int i=0; i<a.size(); ++i)
    {
      it ii=i*N;
      for(int j=0; j<a[i].size(); j++)
        m_data[ii] = a[i][j];
    }
  }

  void add(const std::string& c)
  {
    int num=m_data.size();
    m_data.insert(m_data.end(), N, ' ');
    for(int j=0; j<c.size(); j++,num++)
      m_data[num] = c[j];
  }

  void add(const char* c)
  {
    std::string ctmp(c);
    add(ctmp);
  }

  const char* c_str() const
  {
    return  &(m_data[0]);
  }

};
#endif
