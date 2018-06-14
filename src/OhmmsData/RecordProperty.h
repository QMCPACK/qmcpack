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
    
    


#ifndef OHMMS_RECORDPROPERTYBASE_H__
#define OHMMS_RECORDPROPERTYBASE_H__

#include "OhmmsData/OhmmsElementBase.h"
#include <algorithm>
#include <vector>
#include <string>
#include <fstream>

/** abstract base class to recrod any properties
 */
struct RecordProperty
{

  int stride;
  std::string FileName;
  RecordProperty(): stride(-1),FileName("default") { }

  virtual ~RecordProperty() { }
  virtual void reset(const char* fileroot, bool append=false) = 0;
  virtual void report(int i) = 0;
  virtual void finalize() = 0;
  virtual bool put(xmlNodePtr cur) =0;
};


class RecordPropertyList
{

public:

  typedef std::vector<RecordProperty*> RecordList_t;

  RecordPropertyList() {}
  virtual ~RecordPropertyList()
  {
    clear();
  }

  inline void addRecord(RecordProperty* a)
  {
    Properties.push_back(a);
  }

  inline void setstride(int n)
  {
    RecordList_t::iterator it= Properties.begin();
    while(it != Properties.end())
    {
      (*it)->stride = n;
      it++;
    }
  }

  inline void report(int iter)
  {
    RecordList_t::iterator it= Properties.begin();
    while(it != Properties.end())
    {
      (*it)->report(iter);
      it++;
    }
  }

  inline void finalize()
  {
    RecordList_t::iterator it= Properties.begin();
    while(it != Properties.end())
    {
      (*it)->finalize();
      it++;
    }
  }
  inline void clear()
  {
    std::vector<RecordProperty* >::iterator it = Properties.begin();
    while(it != Properties.end())
    {
      delete (*it);
      it++;
    }
    Properties.clear();
  }

protected:

  std::vector<RecordProperty* > Properties;
};


/** Vectorized record engine for scalar properties
 *
 * A series of values with the name is recorded as a table of multiple columns.
 */
template<class T>
struct RecordNamedProperty: public RecordProperty
{

  std::ostream *OutStream;
  std::vector<T> Values;
  std::vector<std::string> Names;

  RecordNamedProperty(): OutStream(0)
  {
    Values.reserve(20);
    Names.reserve(20);
  }

  explicit RecordNamedProperty(int n)
  {
    Values.resize(n,T());
    Names.resize(n);
  }

  RecordNamedProperty(const RecordNamedProperty<T>& a):
    OutStream(0),Values(a.Values),Names(a.Names)
  { }

  ~RecordNamedProperty()
  {
    if(OutStream)
      delete OutStream;
  }

  void clear()
  {
    Names.clear();
    Values.clear();
  }

  inline T operator[](int i) const
  {
    return Values[i];
  }
  inline T& operator[](int i)
  {
    return Values[i];
  }

  ///iterators to use std algorithms
  inline typename std::vector<T>::iterator begin()
  {
    return Values.begin();
  }
  inline typename std::vector<T>::iterator end()
  {
    return Values.end();
  }
  inline typename std::vector<T>::const_iterator begin() const
  {
    return Values.begin();
  }
  inline typename std::vector<T>::const_iterator end() const
  {
    return Values.end();
  }

  //inline int add(const char* aname) {
  inline int add(const std::string& aname)
  {
    int i=0;
    while(i<Names.size())
    {
      if(Names[i]== aname)
        return i;
      i++;
    }
    Names.push_back(aname);
    Values.push_back(T());
    return Names.size()-1;
  }

  /** add multiple items with the aroot
   * @param aroot root
   * @param first first index
   * @param last last index
   * @return the last valid index
   */
  inline int append(const std::string& aroot, int first, int last)
  {
    for(int i=first; i<last; ++i)
    {
      std::ostringstream o;
      o<<aroot<<i;
      Names.push_back(o.str());
      Values.push_back(T());
    }
    return Values.size();
  }

  inline int size() const
  {
    return Names.size();
  }

  inline void setValues(T v)
  {
    for(int i=0; i<Values.size(); i++)
      Values[i] = v;
  }

  inline void resize(int n)
  {
    std::vector<T> a = Values;
    std::vector<std::string> b = Names;
    Values.resize(n,T());
    for(int i=0; i<a.size(); i++)
      Values[i] = a[i];
    //std::copy_n(a.begin(), a.size(), Values.begin());
    Names.resize(n);
    for(int i=0; i<a.size(); i++)
      Names[i] = b[i];
    //std::copy_n(b.begin(), b.size(), Name.begin());
  }

  ///implement virtual functions
  inline void reset(const char* fileroot, bool append=false)
  {
    if(OutStream)
      delete OutStream;
    if(append)
    {
      OutStream = new std::ofstream(fileroot,std::ios_base::app);
    }
    else
    {
      OutStream = new std::ofstream(fileroot);
    }
    if(!append)
    {
      OutStream->setf(std::ios::left,std::ios::adjustfield);
      *OutStream << "#   ";
      for(int i=0; i<Names.size(); i++)
        (*OutStream) << std::setw(15) << Names[i].c_str();
      (*OutStream) << std::endl;
    }
    OutStream->setf(std::ios::scientific, std::ios::floatfield);
    OutStream->setf(std::ios::right,std::ios::adjustfield);
  }

  inline void report(int iter)
  {
    if(stride && iter%stride == 0)
    {
      for(int i=0; i<Values.size(); i++)
        (*OutStream) << std::setw(15) << Values[i];
      (*OutStream) << std::endl;
    }
  }

  void finalize() {}
  bool put(xmlNodePtr cur);
};

#if defined(HAVE_LIBXML2)
template<class T>
bool RecordNamedProperty<T>::put(xmlNodePtr cur)
{
  xmlAttrPtr att = cur->properties;
  while(att != NULL)
  {
    std::string aname((const char*)(att->name));
    if(aname == "stride")
    {
      stride = atoi((const char*)(att->children->content));
    }
    att = att->next;
  }
  return true;
}
#endif
#endif
