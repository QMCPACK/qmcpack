//////////////////////////////////////////////////////////////////
// (c) Copyright 1998-2002 by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
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
#ifndef OHMMS_FORTRANSTRING_CONVERSION_H
#define OHMMS_FORTRANSTRING_CONVERSION_H

#include <stdio.h>
#include <vector>
#include <string>
#include <iostream>
using std::vector;
using std::string;

/*!\class template<int N> struct FortranString
 *\brief A wrapper class to pass a string to fortran codes
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

  vector<char> m_data;

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

  void set(const vector<string>& a)
  {
    resize(a.size());
    for(int i=0; i<a.size(); ++i)
    {
      it ii=i*N;
      for(int j=0; j<a[i].size(); j++)
        m_data[ii] = a[i][j];
    }
  }

  void add(const string& c)
  {
    int num=m_data.size();
    m_data.insert(m_data.end(), N, ' ');
    for(int j=0; j<c.size(); j++,num++)
      m_data[num] = c[j];
  }

  void add(const char* c)
  {
    string ctmp(c);
    add(ctmp);
  }

  const char* c_str() const
  {
    return  &(m_data[0]);
  }

};
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/
