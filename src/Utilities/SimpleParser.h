//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


#ifndef OHMMS_SIMPLEPARSER_H
#define OHMMS_SIMPLEPARSER_H

#include <cstdlib>
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <list>
#include <sstream>

char* readLine(char *s, int max, std::istream &fp);

int getwords(std::vector<std::string>& slist, std::istream &fp, int dummy = 0, const std::string &extra_tokens="");
int getwords(std::vector<std::string>& slist, std::istream &fp, std::string& aline );
int getwords(std::vector<std::string>& slist,std::istream& fpos, const char* field, const char* terminate);
int getwords(std::vector<std::string>& slist,std::istream& fpos, const char* terminate);
int getXwords(std::vector<std::string>& slist, std::istream &fp);
int getXwords(std::vector<std::string>& slist,std::istream& fpos, const char* terminate);

unsigned parsewords(const char *inbuf, std::vector<std::string>& slist, const std::string &extra_tokens="");
unsigned parsewords(const char *inbuf, std::list<std::string>& slist);

void readXmol( std::istream& ,double*,int);

struct OhmmsAsciiParser
{

  static const int bufferSize=200;
  char dbuffer[bufferSize];
  std::vector<std::string> currentWords;

  inline void skiplines(std::istream& is, int n)
  {
    while(n>0)
    {
      is.getline( dbuffer, bufferSize);
      --n;
    }
  }
  template<class T>
  inline void getValue(std::istream& is, T& aval)
  {
    is.getline( dbuffer,bufferSize);
    std::istringstream a(dbuffer);
    a>> aval;
  }

  template<class T1, class T2>
  inline void getValue(std::istream& is, T1& aval, T2& bval)
  {
    is.getline( dbuffer,bufferSize);
    std::istringstream a(dbuffer);
    a>> aval>>bval;
  }

  template<class IT>
  inline void getValues(std::istream& is, IT first, IT last)
  {
    while(first != last)
    {
      is.getline( dbuffer,bufferSize);
      std::istringstream a(dbuffer);
      while(first != last && a >> *first)
      {
        first++;
      }
    }
  }

  int search(std::istream& is, const std::string& keyword)
  {
    bool notfound = true;
    while(notfound)
    {
      std::string aline;
      getline(is,aline,'\n');
      if(! is)
      {
        std::cout << "KEYWORD " << keyword << " : NOT FOUND. " << std::endl;
        abort();
      }
      if(aline.find(keyword) < aline.size())
      {
        notfound = false;
      }
    }
    return 1;
  }

  int search(std::istream& is, const std::string& keyword, std::string& the_line)
  {
    bool notfound = true;
    while(notfound)
    {
      std::string aline;
      getline(is,aline,'\n');
      if(! is)
      {
        std::cout << "KEYWORD " << keyword << " : NOT FOUND. " << std::endl;
        abort();
      }
      if(aline.find(keyword) < aline.size())
      {
        notfound = false;
        the_line = aline;
      }
    }
    return 1;
  }

  bool lookFor(std::istream& is, const std::string& keyword)
  {
    bool notfound = true;
    while(notfound)
    {
      std::string aline;
      getline(is,aline,'\n');
      if(aline.find(keyword) != std::string::npos)
        // < aline.size()) {
      {
        notfound = false;
      }
      //if(! is){
      if(is.eof())
      {
        return false;
      }
    }
    return true;
  }

  bool lookFor(std::istream& is, const std::string& keyword, std::string& the_line)
  {
    bool notfound = true;
    while(notfound)
    {
      std::string aline;
      getline(is,aline,'\n');
      if(aline.find(keyword) != std::string::npos)
        // < aline.size()) {
      {
        notfound = false;
        the_line = aline;
      }
      //if(! is){
      if(is.eof())
      {
        return false;
      }
    }
    return true;
  }

};

#endif

