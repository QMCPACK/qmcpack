//////////////////////////////////////////////////////////////////
// (c) Copyright 2003  by Jeongnim Kim
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
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef OHMMS_RECORDPROPERTYBASE_H__
#define OHMMS_RECORDPROPERTYBASE_H__

#include "OhmmsData/OhmmsElementBase.h"
#include <algorithm>
#include <vector>
#include <string>
#include <fstream>

struct RecordProperty {
  
  int stride;
  std::string FileName;
  RecordProperty(): stride(-1),FileName("default"){ }
  
  virtual ~RecordProperty() { }
  virtual void reset(const char* fileroot) = 0;
  virtual void report(int i) = 0;
  virtual void finalize() = 0;
  virtual bool put(xmlNodePtr cur) =0;
};


class RecordPropertyList {

 public:

  typedef std::vector<RecordProperty*> RecordList_t;
  
  RecordPropertyList(){}
  virtual ~RecordPropertyList(){ clear();}
  
  inline void addRecord(RecordProperty* a) {
    Properties.push_back(a);
  }
  
  inline void setstride(int n) {
    RecordList_t::iterator it= Properties.begin();
    while(it != Properties.end()) {
      (*it)->stride = n; it++;
    }
  }
  
  inline void report(int iter) {
    RecordList_t::iterator it= Properties.begin();
    while(it != Properties.end()) {
      (*it)->report(iter); it++;
    }
  }
  
  inline void finalize() {
    RecordList_t::iterator it= Properties.begin();
    while(it != Properties.end()) {
      (*it)->finalize(); it++;
    }
  }
  inline void clear() {
    std::vector<RecordProperty* >::iterator it = Properties.begin();
    while(it != Properties.end()) {delete (*it); it++;}
    Properties.clear();
  }
  
 protected:
  
  std::vector<RecordProperty* > Properties;
};


template<class T>
struct RecordNamedProperty: public RecordProperty {

  std::ostream *OutStream;
  std::vector<T> Values;
  std::vector<std::string> Name;

  RecordNamedProperty(): OutStream(NULL) {
    Values.reserve(20); 
    Name.reserve(20);
  }

  explicit RecordNamedProperty(int n) { 
    Values.resize(n,T());
    Name.resize(n);
  }

  ~RecordNamedProperty() { 
    if(OutStream) delete OutStream;
  }

  inline T operator[](int i) const { return Values[i];}
  inline T& operator[](int i) { return Values[i];}

  inline int add(const char* aname) {
    int i=0;
    while(i<Name.size()) {
      if(Name[i]== aname) return i;
      i++;
    }
    Name.push_back(aname);
    Values.push_back(T());
    return Name.size()-1;
  }

  inline int size() const { return Name.size();}

  inline void resize(int n) {
    std::vector<T> a = Values;
    std::vector<std::string> b = Name;
    Values.resize(n,T()); 
    for(int i=0; i<a.size(); i++) Values[i] = a[i];
    //std::copy_n(a.begin(), a.size(), Values.begin());
    Name.resize(n);  
    for(int i=0; i<a.size(); i++) Name[i] = b[i];
    //std::copy_n(b.begin(), b.size(), Name.begin());
  }


  ///implement virtual functions
  inline void reset(const char* fileroot) { 
    if(OutStream) delete OutStream;
    OutStream = new std::ofstream(fileroot);
    OutStream->setf(std::ios::scientific, std::ios::floatfield);
    OutStream->setf(std::ios::left,std::ios::adjustfield);
    *OutStream << "#   ";
    for(int i=0; i<Name.size(); i++) 
      (*OutStream) << setw(15) << Name[i].c_str();
    (*OutStream) << endl;
    OutStream->setf(std::ios::right,std::ios::adjustfield);
  }

  inline void report(int iter){ 
    if(stride && iter%stride == 0) {
      for(int i=0; i<Values.size(); i++) 
	(*OutStream) << setw(15) << Values[i];
      (*OutStream) << endl;
    }
  }

  void finalize() {}
  bool put(xmlNodePtr cur);
};

#if defined(HAVE_LIBXML2)
template<class T>
bool RecordNamedProperty<T>::put(xmlNodePtr cur){

  xmlAttrPtr att = cur->properties;
  while(att != NULL) {
    std::string aname((const char*)(att->name));
    if(aname == "stride") {
      stride = atoi((const char*)(att->children->content));
    }
    att = att->next;
  }
  return true;
}
#endif
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
