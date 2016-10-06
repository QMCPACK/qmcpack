//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
    
    




#ifndef STRING_UTILS_H
#define STRING_UTILS_H


#include<cstdio>
#include<Configuration.h>

namespace qmcplusplus
{


inline std::string strip(const std::string& s)
{
  int start=s.length();
  int end=0;
  int i;
  for(i=0; i<s.length(); i++)
  {
    if(s[i]!=' '&&s[i]!='\n'&&s[i]!='\t')
    {
      start=i;
      break;
    }
  }
  for(i=s.length()-1; i>0; i--)
  {
    if(s[i]!=' '&&s[i]!='\n'&&s[i]!='\t')
    {
      end=i;
      break;
    }
  }
  //app_log()<<"strip got '"<<s<<"'"<< std::endl;
  //app_log()<<"start,end "<<start<<","<<end<<" "<<s[start]<<" "<<s[end]<< std::endl;
  //app_log()<<"returning '"<<s.substr(start,end-start+1)<<"'"<< std::endl;
  return s.substr(start,end-start+1);
}


inline bool whitespace(char c)
{
  return (c==' ' || c=='\n' || c=='\t');
}


inline std::vector<std::string> split(const std::string& s)
{
  std::vector<std::string> tokens;
  int i=0;
  while(i<s.length())
  {
    while(i<s.length() && whitespace(s[i]))
      i++;
    int start = i;
    while(i<s.length() && !whitespace(s[i]))
      i++;
    int end = i;
    int len = end-start;
    if(len>0)
      tokens.push_back(s.substr(start,len));
  }
  return tokens;
}


inline std::vector<std::string> split(const std::string& s, const std::string& pattern)
{
  int sloc=0;
  int eloc;
  int plen=pattern.length();
  std::string ss;
  std::vector<std::string> tokens;
  //app_log() << "split got string:" << std::endl<<"'"<<s<<"'"<< std::endl;
  while(true)
  {
    eloc=s.find(pattern,sloc);
    if(eloc!=std::string::npos)
    {
      ss=s.substr(sloc,eloc-sloc);
      if(ss!="")
      {
        //app_log()<<"  adding token: "<< std::endl;
        //app_log()<<"    '"<< ss <<"'" << std::endl;
        tokens.push_back(ss);
      }
      sloc=eloc+plen;
    }
    else
    {
      eloc=s.length();
      ss=s.substr(sloc,eloc-sloc);
      if(ss!="")
      {
        //app_log()<<"  adding token: "<< std::endl;
        //app_log()<<"    '"<< ss <<"'" << std::endl;
        tokens.push_back(ss);
      }
      break;
    }
  }
  return tokens;
}

inline int string2int(const std::string& s)
{
  return atoi(s.c_str());
}

inline double string2real(const std::string& s)
{
  return atof(s.c_str());
}

inline std::string int2string(const int& i)
{
  std::stringstream ss;
  ss<<i;
  return ss.str();
}

inline std::string real2string(const double& r)
{
  std::stringstream ss;
  ss<<r;
  return ss.str();
}

inline bool string2bool(const std::string& s)
{
  if(s=="true" || s=="yes" || s=="1")
  {
    return true;
  }
  else
    if(s=="false" || s=="no" || s=="0")
    {
      return false;
    }
    else
    {
      APP_ABORT("string2bool received non-boolean string: "+s);
      return false;
    }
}


//strings for input (OhmmsAttributeSet)
struct astring
{
  std::string s;
};
inline std::istream& operator>>( std::istream& is,astring& rhs)
{
  char buf[256];
  is.getline(buf,256);
  rhs.s.assign(buf);
  return is;
}
inline std::ostream& operator<<(std::ostream& os, const astring& rhs)
{
  os<<rhs.s<< std::endl;
  return os;
}

}

#endif
