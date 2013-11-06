
#ifndef STRING_UTILS_H
#define STRING_UTILS_H


#include<cstdio>
#include<Configuration.h>

namespace qmcplusplus
{

using namespace std;

inline string strip(const string& s)
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
  //app_log()<<"strip got '"<<s<<"'"<<endl;
  //app_log()<<"start,end "<<start<<","<<end<<" "<<s[start]<<" "<<s[end]<<endl;
  //app_log()<<"returning '"<<s.substr(start,end-start+1)<<"'"<<endl;
  return s.substr(start,end-start+1);
}


inline bool whitespace(char c)
{
  return (c==' ' || c=='\n' || c=='\t');
}


inline vector<string> split(const string& s)
{
  vector<string> tokens;
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


inline vector<string> split(const string& s, const string& pattern)
{
  int sloc=0;
  int eloc;
  int plen=pattern.length();
  string ss;
  vector<string> tokens;
  //app_log() << "split got string:" <<endl<<"'"<<s<<"'"<<endl;
  while(true)
  {
    eloc=s.find(pattern,sloc);
    if(eloc!=string::npos)
    {
      ss=s.substr(sloc,eloc-sloc);
      if(ss!="")
      {
        //app_log()<<"  adding token: "<<endl;
        //app_log()<<"    '"<< ss <<"'" <<endl;
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
        //app_log()<<"  adding token: "<<endl;
        //app_log()<<"    '"<< ss <<"'" <<endl;
        tokens.push_back(ss);
      }
      break;
    }
  }
  return tokens;
}

inline int string2int(const string& s)
{
  return atoi(s.c_str());
}

inline double string2real(const string& s)
{
  return atof(s.c_str());
}

inline string int2string(const int& i)
{
  stringstream ss;
  ss<<i;
  return ss.str();
}

inline string real2string(const double& r)
{
  stringstream ss;
  ss<<r;
  return ss.str();
}

inline bool string2bool(const string& s)
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
  string s;
};
inline istream& operator>>(istream& is,astring& rhs)
{
  char buf[256];
  is.getline(buf,256);
  rhs.s.assign(buf);
  return is;
}
inline ostream& operator<<(ostream& os, const astring& rhs)
{
  os<<rhs.s<<endl;
  return os;
}

}

#endif
